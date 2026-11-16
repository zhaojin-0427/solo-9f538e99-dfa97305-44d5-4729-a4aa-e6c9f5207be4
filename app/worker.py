"""后台计算 worker：认领任务 → 计算 → 单事务写入报告。

并发语义：
- 报告与任务绑定固定的 (fasta_version, csv_version)，是否"当前"由版本派生，
  因此迟到完成的结果永远不会覆盖更新版本对应的结论；
- 报告 INSERT 与任务状态 UPDATE 在同一事务提交，失败不留半份报告。
"""
from __future__ import annotations

import json
import time
import traceback

from .db import connect
from .normalize import ValidationError, compute_report


def claim_job(db_path: str):
    """原子认领最早的 pending 任务；没有则返回 None。"""
    conn = connect(db_path)
    try:
        row = conn.execute(
            "SELECT id FROM jobs WHERE status='pending' ORDER BY id LIMIT 1"
        ).fetchone()
        if row is None:
            return None
        cur = conn.execute(
            "UPDATE jobs SET status='running', started_at=datetime('now')"
            " WHERE id=? AND status='pending'",
            (row["id"],),
        )
        conn.commit()
        if cur.rowcount != 1:  # 被其他 worker 抢先
            return None
        job = conn.execute("SELECT * FROM jobs WHERE id=?", (row["id"],)).fetchone()
        return dict(job)
    finally:
        conn.close()


def process_job(db_path: str, job: dict, delay: float = 0.0) -> None:
    """执行一个任务并把报告与任务状态在同一事务落库。"""
    conn = connect(db_path)
    try:
        fasta = conn.execute(
            "SELECT content FROM inputs WHERE batch_id=? AND kind='fasta' AND version=?",
            (job["batch_id"], job["fasta_version"]),
        ).fetchone()
        csv_row = conn.execute(
            "SELECT content FROM inputs WHERE batch_id=? AND kind='csv' AND version=?",
            (job["batch_id"], job["csv_version"]),
        ).fetchone()

        if delay > 0:
            time.sleep(delay)  # 模拟耗时计算；此期间更正照常提交

        try:
            if fasta is None or csv_row is None:
                raise ValidationError("任务对应的输入版本不存在")
            result = compute_report(fasta["content"], csv_row["content"])
            status, error, payload = "ok", None, json.dumps(result, ensure_ascii=False)
            chrom, ref_len = result["chrom"], result["ref_len"]
        except ValidationError as exc:
            status, error, payload, chrom, ref_len = "failed", str(exc), None, None, None
        except Exception as exc:  # 非预期错误同样整批失败，不外泄堆栈
            status, error = "failed", f"内部错误: {exc}"
            payload, chrom, ref_len = None, None, None

        job_status = "done" if status == "ok" else "failed"
        with conn:  # 单事务：报告 + 任务状态一起提交，失败整体回滚
            conn.execute(
                "INSERT INTO reports(job_id, batch_id, fasta_version, csv_version,"
                " status, error, chrom, ref_len, payload, created_at)"
                " VALUES (?,?,?,?,?,?,?,?,?, datetime('now'))",
                (job["id"], job["batch_id"], job["fasta_version"], job["csv_version"],
                 status, error, chrom, ref_len, payload),
            )
            conn.execute(
                "UPDATE jobs SET status=?, finished_at=datetime('now'), error=? WHERE id=?",
                (job_status, error, job["id"]),
            )
    finally:
        conn.close()


def run_pending_jobs(db_path: str, delay: float = 0.0, limit=None) -> int:
    """依次处理 pending 任务，返回处理数（测试与 worker 循环共用）。"""
    n = 0
    while limit is None or n < limit:
        job = claim_job(db_path)
        if job is None:
            break
        process_job(db_path, job, delay)
        n += 1
    return n


def worker_loop(db_path: str, delay: float, poll: float) -> None:
    while True:
        try:
            if run_pending_jobs(db_path, delay, limit=1) == 0:
                time.sleep(poll)
        except Exception:
            traceback.print_exc()
            time.sleep(poll)
