"""HTTP 路由：批次创建、输入更正、报告查看。"""
from __future__ import annotations

import json

from flask import (Blueprint, Response, abort, flash, redirect,
                   render_template, request, url_for)

from .db import get_db

bp = Blueprint("main", __name__)

KIND_LABEL = {"fasta": "FASTA", "csv": "CSV"}


# ---------------------------------------------------------------- 辅助

def _batch_or_404(db, bid):
    batch = db.execute("SELECT * FROM batches WHERE id=?", (bid,)).fetchone()
    if batch is None:
        abort(404)
    return batch


def _current_versions(db, bid):
    rows = db.execute(
        "SELECT kind, MAX(version) AS v FROM inputs WHERE batch_id=? GROUP BY kind",
        (bid,),
    ).fetchall()
    d = {r["kind"]: r["v"] for r in rows}
    return d.get("fasta"), d.get("csv")


def _enqueue(db, bid):
    """按当前最新版本组合排队一个计算任务（相同版本组合的未完成任务不重复排）。"""
    fv, cv = _current_versions(db, bid)
    dup = db.execute(
        "SELECT 1 FROM jobs WHERE batch_id=? AND fasta_version=? AND csv_version=?"
        " AND status IN ('pending','running')",
        (bid, fv, cv),
    ).fetchone()
    if dup is None:
        db.execute(
            "INSERT INTO jobs(batch_id, fasta_version, csv_version, status, created_at)"
            " VALUES (?,?,?,'pending', datetime('now'))",
            (bid, fv, cv),
        )


def _read_upload(field):
    """读取上传文件并解码为文本；失败时 flash 并返回 None。"""
    f = request.files.get(field)
    if f is None or not f.filename:
        flash("请选择要上传的文件")
        return None
    try:
        return f.read().decode("utf-8-sig")
    except UnicodeDecodeError:
        flash("文件必须是 UTF-8 编码的文本")
        return None


# ---------------------------------------------------------------- 页面

@bp.get("/")
def index():
    db = get_db()
    batches = db.execute("SELECT * FROM batches ORDER BY id DESC").fetchall()
    return render_template("index.html", batches=batches)


@bp.post("/batches")
def create_batch():
    name = request.form.get("name", "").strip() or "未命名批次"
    fasta_text = _read_upload("fasta")
    csv_text = _read_upload("csv")
    if fasta_text is None or csv_text is None:
        return redirect(url_for("main.index"))

    db = get_db()
    cur = db.execute(
        "INSERT INTO batches(name, created_at) VALUES (?, datetime('now'))", (name,)
    )
    bid = cur.lastrowid
    for kind, text in (("fasta", fasta_text), ("csv", csv_text)):
        db.execute(
            "INSERT INTO inputs(batch_id, kind, version, content, created_at)"
            " VALUES (?,?,1,?, datetime('now'))",
            (bid, kind, text),
        )
    _enqueue(db, bid)
    db.commit()
    flash(f"批次 #{bid} 已创建（FASTA v1 / CSV v1），计算已排队")
    return redirect(url_for("main.batch_detail", bid=bid))


@bp.post("/batches/<int:bid>/inputs/<kind>")
def correct_input(bid, kind):
    """整份更正某一份输入；计算进行中也可提交。"""
    if kind not in KIND_LABEL:
        abort(404)
    db = get_db()
    _batch_or_404(db, bid)
    text = _read_upload("file")
    if text is None:
        return redirect(url_for("main.batch_detail", bid=bid))

    cur = db.execute(
        "SELECT MAX(version) AS v FROM inputs WHERE batch_id=? AND kind=?",
        (bid, kind),
    ).fetchone()["v"]
    db.execute(
        "INSERT INTO inputs(batch_id, kind, version, content, created_at)"
        " VALUES (?,?,?,?, datetime('now'))",
        (bid, kind, cur + 1, text),
    )
    _enqueue(db, bid)
    db.commit()
    flash(f"{KIND_LABEL[kind]} 第 {cur + 1} 版已提交，已按最新版本组合排队重新计算")
    return redirect(url_for("main.batch_detail", bid=bid))


@bp.get("/batches/<int:bid>")
def batch_detail(bid):
    db = get_db()
    batch = _batch_or_404(db, bid)
    fv, cv = _current_versions(db, bid)
    report = db.execute(
        "SELECT * FROM reports WHERE batch_id=? AND fasta_version=? AND csv_version=?"
        " ORDER BY id DESC LIMIT 1",
        (bid, fv, cv),
    ).fetchone()
    computing = db.execute(
        "SELECT COUNT(*) AS n FROM jobs WHERE batch_id=? AND fasta_version=?"
        " AND csv_version=? AND status IN ('pending','running')",
        (bid, fv, cv),
    ).fetchone()["n"] > 0
    reports = db.execute(
        "SELECT * FROM reports WHERE batch_id=? ORDER BY id DESC LIMIT 50", (bid,)
    ).fetchall()
    jobs = db.execute(
        "SELECT * FROM jobs WHERE batch_id=? ORDER BY id DESC LIMIT 20", (bid,)
    ).fetchall()
    payload = (json.loads(report["payload"])
               if report and report["status"] == "ok" else None)
    return render_template(
        "batch.html", batch=batch, fv=fv, cv=cv, report=report, payload=payload,
        computing=computing, reports=reports, jobs=jobs,
    )


@bp.get("/batches/<int:bid>/reports/<int:rid>")
def report_detail(bid, rid):
    db = get_db()
    batch = _batch_or_404(db, bid)
    report = db.execute(
        "SELECT * FROM reports WHERE id=? AND batch_id=?", (rid, bid)
    ).fetchone()
    if report is None:
        abort(404)
    fv, cv = _current_versions(db, bid)
    is_current = (report["fasta_version"] == fv and report["csv_version"] == cv)
    payload = (json.loads(report["payload"])
               if report["status"] == "ok" else None)
    return render_template(
        "report.html", batch=batch, report=report, payload=payload,
        is_current=is_current, fv=fv, cv=cv,
    )


@bp.get("/batches/<int:bid>/reports/<int:rid>/vcf")
def report_vcf(bid, rid):
    db = get_db()
    report = db.execute(
        "SELECT * FROM reports WHERE id=? AND batch_id=?", (rid, bid)
    ).fetchone()
    if report is None or report["status"] != "ok":
        abort(404)
    vcf = json.loads(report["payload"])["vcf"]
    return Response(
        vcf,
        content_type="text/plain; charset=utf-8",
        headers={"Content-Disposition":
                 f"attachment; filename=batch{bid}_report{rid}.vcf"},
    )
