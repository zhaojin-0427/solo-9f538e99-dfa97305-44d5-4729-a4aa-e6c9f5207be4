"""应用层测试：批次流程、版本化更正、迟到结果不覆盖、失败原子性。"""
import io
import json
import time

import pytest

from app import create_app
from app.db import connect
from app.worker import claim_job, process_job, run_pending_jobs

FASTA = ">chr1\n" + "A" * 20 + "\n"
CSV_OK = "ID,POS,REF,ALT\ni1,3,A,AA\ni2,9,A,AA\ns1,5,A,C\n"
CSV_BAD = "ID,POS,REF,ALT\nbad,5,C,G\n"          # REF 与参考不符
CSV_OK2 = "ID,POS,REF,ALT\nd1,10,AA,A\n"         # 更正后的 CSV


@pytest.fixture()
def app(tmp_path):
    return create_app({
        "TESTING": True,
        "DATABASE": str(tmp_path / "t.db"),
        "WORKER_ENABLED": False,     # 测试手动驱动 worker
        "COMPUTE_DELAY_SECONDS": 0,
    })


@pytest.fixture()
def client(app):
    return app.test_client()


def _upload(client, fasta=FASTA, csv_text=CSV_OK, name="t"):
    return client.post("/batches", data={
        "name": name,
        "fasta": (io.BytesIO(fasta.encode()), "ref.fa"),
        "csv": (io.BytesIO(csv_text.encode()), "vars.csv"),
    }, follow_redirects=False)


def _correct(client, bid, kind, text):
    return client.post(f"/batches/{bid}/inputs/{kind}", data={
        "file": (io.BytesIO(text.encode()), f"new.{kind}"),
    })


def _db(app):
    return connect(app.config["DATABASE"])


def _reports(app, bid=1):
    conn = _db(app)
    rows = conn.execute(
        "SELECT * FROM reports WHERE batch_id=? ORDER BY id", (bid,)).fetchall()
    conn.close()
    return rows


# ---------------------------------------------------------------- 基本流程

def test_create_compute_and_current_report(app, client):
    resp = _upload(client)
    assert resp.status_code == 302
    assert run_pending_jobs(app.config["DATABASE"]) == 1

    reports = _reports(app)
    assert len(reports) == 1 and reports[0]["status"] == "ok"
    payload = json.loads(reports[0]["payload"])
    assert len(payload["groups"]) == 2            # 两个插入归并 + 一个 SNV
    ins = [g for g in payload["groups"] if g["kind"] == "INS"][0]
    assert ins["ids"] == ["i1", "i2"] and ins["pos"] == 1

    page = client.get("/batches/1").data.decode()
    assert "当前报告" in page and "i1, i2" in page


def test_invalid_row_fails_whole_batch_without_partial_report(app, client):
    _upload(client, csv_text=CSV_BAD)
    run_pending_jobs(app.config["DATABASE"])
    reports = _reports(app)
    assert len(reports) == 1 and reports[0]["status"] == "failed"
    assert reports[0]["payload"] is None          # 失败不留半份报告
    assert "与参考序列不符" in reports[0]["error"]
    page = client.get("/batches/1").data.decode()
    assert "失败" in page and "与参考序列不符" in page


def test_vcf_download(app, client):
    _upload(client)
    run_pending_jobs(app.config["DATABASE"])
    resp = client.get("/batches/1/reports/1/vcf")
    assert resp.status_code == 200
    body = resp.data.decode()
    assert body.startswith("##fileformat=VCFv4.2")
    assert "chr1\t1\ti1;i2\tA\tAA" in body


# ---------------------------------------------------------------- 版本与更正

def test_correction_creates_new_version_and_supersedes_old_report(app, client):
    _upload(client)
    run_pending_jobs(app.config["DATABASE"])

    _correct(client, 1, "csv", CSV_OK2)
    run_pending_jobs(app.config["DATABASE"])

    reports = _reports(app)
    assert len(reports) == 2
    old, new = reports
    assert (old["fasta_version"], old["csv_version"]) == (1, 1)
    assert (new["fasta_version"], new["csv_version"]) == (1, 2)

    # 旧报告可查，但不再是当前结论
    old_page = client.get("/batches/1/reports/1").data.decode()
    assert "已过期" in old_page
    new_page = client.get("/batches/1/reports/2").data.decode()
    assert "当前" in new_page
    batch_page = client.get("/batches/1").data.decode()
    assert "CSV v2" in batch_page and "d1" in batch_page


def test_fasta_correction_keeps_csv_version(app, client):
    _upload(client)
    run_pending_jobs(app.config["DATABASE"])
    _correct(client, 1, "fasta", ">chr1\n" + "A" * 30 + "\n")
    run_pending_jobs(app.config["DATABASE"])
    reports = _reports(app)
    assert (reports[-1]["fasta_version"], reports[-1]["csv_version"]) == (2, 1)


def test_late_result_cannot_overwrite_newer_version(app, client):
    """任务计算期间提交更正：迟到完成的旧版本报告不能成为当前结论。"""
    _upload(client)
    db_path = app.config["DATABASE"]

    job1 = claim_job(db_path)                     # 任务 1 开始计算（旧版本）
    assert job1 is not None
    _correct(client, 1, "csv", CSV_OK2)           # 计算期间更正 CSV → v2 + 新任务
    process_job(db_path, job1)                    # 旧任务迟到完成

    conn = _db(app)
    cur = conn.execute(
        "SELECT MAX(version) v FROM inputs WHERE batch_id=1 AND kind='csv'").fetchone()
    assert cur["v"] == 2
    # 当前版本 (1,2) 尚无成功报告 → 批次页显示计算中，而非迟到的旧结论
    page = client.get("/batches/1").data.decode()
    assert "计算中" in page
    conn.close()

    run_pending_jobs(db_path)                     # 新版本任务完成
    reports = _reports(app)
    assert len(reports) == 2
    page = client.get("/batches/1").data.decode()
    assert "d1" in page                           # 当前结论来自新版本
    old_page = client.get("/batches/1/reports/1").data.decode()
    assert "已过期" in old_page


def test_failed_correction_leaves_no_partial_and_can_be_fixed(app, client):
    _upload(client)
    run_pending_jobs(app.config["DATABASE"])
    _correct(client, 1, "csv", CSV_BAD)           # 更正为非法 CSV
    run_pending_jobs(app.config["DATABASE"])
    reports = _reports(app)
    assert reports[-1]["status"] == "failed" and reports[-1]["payload"] is None
    page = client.get("/batches/1").data.decode()
    assert "失败" in page

    _correct(client, 1, "csv", CSV_OK2)           # 再次更正恢复正常
    run_pending_jobs(app.config["DATABASE"])
    reports = _reports(app)
    assert reports[-1]["status"] == "ok"
    assert (reports[-1]["fasta_version"], reports[-1]["csv_version"]) == (1, 3)


def test_old_reports_remain_viewable(app, client):
    _upload(client)
    run_pending_jobs(app.config["DATABASE"])
    _correct(client, 1, "csv", CSV_OK2)
    run_pending_jobs(app.config["DATABASE"])
    resp = client.get("/batches/1/reports/1")
    assert resp.status_code == 200 and "i1" in resp.data.decode()


# ---------------------------------------------------------------- 并发集成（真实 worker 线程）

def test_correction_accepted_while_computing(tmp_path):
    db_path = str(tmp_path / "live.db")
    app = create_app({
        "TESTING": True,
        "DATABASE": db_path,
        "WORKER_ENABLED": True,
        "COMPUTE_DELAY_SECONDS": 1.0,   # 模拟耗时计算
        "WORKER_POLL_SECONDS": 0.05,
    })
    client = app.test_client()
    _upload(client)
    time.sleep(0.3)                               # 任务应已进入 running
    _correct(client, 1, "csv", CSV_OK2)           # 计算期间更正仍可提交

    deadline = time.time() + 10
    while time.time() < deadline:
        conn = connect(db_path)
        n = conn.execute("SELECT COUNT(*) c FROM reports").fetchone()["c"]
        conn.close()
        if n >= 2:
            break
        time.sleep(0.1)
    reports = _reports(app)
    assert len(reports) == 2
    assert reports[0]["csv_version"] == 1 and reports[1]["csv_version"] == 2
    page = client.get("/batches/1").data.decode()
    assert "d1" in page                           # 最终当前结论属于新版本
