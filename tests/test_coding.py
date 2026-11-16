"""编码影响分析测试：CDS 拼接/翻译、位点判定、整次拒绝、版本语义。"""
import io
import json

import pytest

from app import create_app
from app.coding import (ValidationError, analyze_coding, first_changed_position,
                        reverse_complement, translate, validate_exons)
from app.db import connect
from app.normalize import compute_report
from app.worker import run_pending_jobs

# ATG | GAA | TTT | TAA，12 bp：M E F *
CDS_REF = "ATGGAATTTTAA"
FASTA = ">chr1\n" + CDS_REF + "\n"


def csv_of(*rows):
    return "ID,POS,REF,ALT\n" + "\n".join(",".join(map(str, r)) for r in rows) + "\n"


def _report(fasta, csv_text):
    return compute_report(fasta, csv_text)


def _site(csv_text, exons, strand="+", fasta=FASTA):
    rep = _report(fasta, csv_text)
    res = analyze_coding(rep["ref_seq"], rep["groups"], exons, strand)
    return res, res["sites"]


# ---------------------------------------------------------------- 基础工具

def test_translate_basic_and_stop():
    assert translate(CDS_REF) == "MEF*"
    assert translate("ATGAAA") == "MK"
    assert translate("ATGNNN") == "MX"


def test_reverse_complement():
    assert reverse_complement("ATGG") == "CCAT"


def test_first_changed_position():
    assert first_changed_position("MEF*", "MKF*") == 2
    assert first_changed_position("MEF*", "MEF*") is None
    assert first_changed_position("MEF*", "MEFK*") == 4


# ---------------------------------------------------------------- 外显子校验（整次拒绝）

def test_exon_out_of_bounds_rejected():
    with pytest.raises(ValidationError, match="越界"):
        validate_exons(12, [{"start": 1, "end": 13}], "+")


def test_exon_start_after_end_rejected():
    with pytest.raises(ValidationError, match="大于终点"):
        validate_exons(12, [{"start": 8, "end": 4}], "+")


def test_exons_overlap_rejected():
    with pytest.raises(ValidationError, match="重叠"):
        validate_exons(12, [{"start": 1, "end": 8}, {"start": 8, "end": 12}], "+")


def test_spliced_length_not_multiple_of_three_rejected():
    with pytest.raises(ValidationError, match="3 的倍数"):
        validate_exons(12, [{"start": 1, "end": 4}, {"start": 5, "end": 8}], "+")


def test_bad_strand_rejected():
    with pytest.raises(ValidationError, match="链方向"):
        validate_exons(12, [{"start": 1, "end": 12}], "x")


def test_more_than_two_exons_rejected():
    with pytest.raises(ValidationError, match="1 至 2 段"):
        validate_exons(12, [{"start": 1, "end": 3}] * 3, "+")


def test_non_numeric_coordinate_rejected():
    with pytest.raises(ValidationError, match="整数"):
        validate_exons(12, [{"start": "", "end": "4"}], "+")


def test_exon_order_is_normalized():
    exons = validate_exons(12, [{"start": 7, "end": 9}, {"start": 1, "end": 3}], "+")
    assert exons == [{"start": 1, "end": 3}, {"start": 7, "end": 9}]


# ---------------------------------------------------------------- 正链位点

def test_snv_missense_on_plus_strand():
    res, sites = _site(csv_of(("s1", 4, "G", "A")), [{"start": 1, "end": 12}])
    s = sites[0]
    assert s["status"] == "coding"
    assert s["ref_protein"] == "MEF*" and s["alt_protein"] == "MKF*"
    assert s["first_change"] == 2
    assert not s["frameshift"] and not s["premature_stop"]


def test_synonymous_snv_reports_no_change():
    # GAA（E）→ GAG（E），第 6 位 A→G
    _, sites = _site(csv_of(("s1", 6, "A", "G")), [{"start": 1, "end": 12}])
    s = sites[0]
    assert s["status"] == "coding"
    assert s["alt_protein"] == "MEF*"
    assert s["first_change"] is None


def test_snv_creating_stop_is_premature():
    # 第 4 位 G→T：密码子 GAA→TAA（终止）
    _, sites = _site(csv_of(("s1", 4, "G", "T")), [{"start": 1, "end": 12}])
    s = sites[0]
    assert s["alt_protein"] == "M*"
    assert s["premature_stop"] is True
    assert s["first_change"] == 2
    assert s["frameshift"] is False


def test_snv_in_intron_is_non_coding():
    fasta = ">chr1\n" + "A" * 6 + CDS_REF + "\n"
    res, sites = _site(csv_of(("s1", 3, "A", "C")),
                       [{"start": 7, "end": 18}], fasta=fasta)
    s = sites[0]
    assert s["status"] == "non_coding"
    assert s["ref_protein"] is None and s["alt_protein"] is None
    assert "不在任何编码外显子" in s["reason"]


def test_insertion_frameshift_and_premature_stop():
    # 在密码子 2 内插入 1 个碱基 → 移码
    _, sites = _site(csv_of(("i1", 4, "G", "GT")), [{"start": 1, "end": 12}])
    s = sites[0]
    assert s["status"] == "coding"
    assert s["frameshift"] is True
    assert s["length_delta"] == 1
    assert s["alt_protein"] == translate(CDS_REF[:4] + "T" + CDS_REF[4:])


def test_insertion_three_bases_no_frameshift():
    _, sites = _site(csv_of(("i1", 4, "G", "GTTT")), [{"start": 1, "end": 12}])
    s = sites[0]
    assert s["status"] == "coding"
    assert s["frameshift"] is False
    assert s["length_delta"] == 3


def test_deletion_of_one_base_is_frameshift():
    # 锚定第 4 位 G，删除第 5 位 A
    _, sites = _site(csv_of(("d1", 4, "GA", "G")), [{"start": 1, "end": 12}])
    assert sites[0]["length_delta"] == -1
    assert sites[0]["frameshift"] is True


def test_deletion_of_three_bases_is_in_frame():
    # 锚定第 3 位 G，删除第 4-6 位 GAA
    _, sites = _site(csv_of(("d1", 3, "GGAA", "G")), [{"start": 1, "end": 12}])
    s = sites[0]
    assert s["status"] == "coding"
    assert s["length_delta"] == -3
    assert s["frameshift"] is False
    assert s["alt_protein"] == "MF*"


EXON_FASTA = ">chr1\n" + CDS_REF[:6] + "TTT" + CDS_REF[6:] + "\n"  # 15 bp，含 3 bp 内含子
EXON_DEF = [{"start": 1, "end": 6}, {"start": 10, "end": 15}]


def test_insertion_at_exon_boundary_is_undecidable():
    # 插入点在第 6/7 位之间：第 7 位属于内含子 → 边界
    res, sites = _site(csv_of(("i1", 6, "A", "AT")), EXON_DEF, fasta=EXON_FASTA)
    assert res["cds_len"] == 12
    s = sites[0]
    assert s["status"] == "boundary"
    assert s["alt_protein"] is None


def test_deletion_spanning_exon_junction_is_undecidable():
    # 删除第 6（外显子）与第 7（内含子）位
    res, sites = _site(csv_of(("d1", 5, "AAT", "A")), EXON_DEF, fasta=EXON_FASTA)
    s = sites[0]
    assert s["status"] == "boundary"
    assert "跨外显子" in s["reason"]


def test_deletion_entirely_intronic_is_non_coding():
    # 删除第 7-8 位（均在内含子）
    res, sites = _site(csv_of(("d1", 6, "ATT", "A")), EXON_DEF, fasta=EXON_FASTA)
    assert sites[0]["status"] == "non_coding"


def test_insertion_at_reference_end_is_boundary():
    # 在外显子末端（也是参考末端）后插入，没有右侧参考碱基可判定
    _, sites = _site(csv_of(("i1", 12, "A", "AT")), [{"start": 1, "end": 12}])
    assert sites[0]["status"] == "boundary"
    assert sites[0]["frameshift"] is False


# ---------------------------------------------------------------- 负链

def test_minus_strand_snv_matches_plus_strand_result():
    plus_res, plus_sites = _site(csv_of(("s1", 4, "G", "A")),
                                [{"start": 1, "end": 12}], "+")
    # 负链使用反向互补参考序列时应得到同一蛋白与同一判定
    rc_fasta = ">chr1\n" + reverse_complement(CDS_REF) + "\n"
    rep = _report(rc_fasta, csv_of(("s1", 9, "C", "T")))  # 对应正链第 4 位 G→A
    res = analyze_coding(rep["ref_seq"], rep["groups"],
                         [{"start": 1, "end": 12}], "-")
    assert res["cds_seq"] == CDS_REF
    assert res["sites"][0]["status"] == "coding"
    assert res["sites"][0]["alt_protein"] == plus_sites[0]["alt_protein"]


def test_minus_strand_insertion_matches_plus_strand_result():
    """同一插入在正链与负链（反向互补参考）上分析结果一致。"""
    _, plus_sites = _site(csv_of(("i1", 4, "G", "GT")),
                          [{"start": 1, "end": 12}], "+")
    rc = reverse_complement(CDS_REF)
    fasta = ">chr1\n" + rc + "\n"
    # 正链第 4/5 位之间插入 T ⇔ 负链第 8/9 位之间插入 A（锚定第 8 位 T）
    rep = _report(fasta, csv_of(("i1", 8, "T", "TA")))
    res = analyze_coding(rep["ref_seq"], rep["groups"],
                         [{"start": 1, "end": 12}], "-")
    s = res["sites"][0]
    assert s["status"] == "coding"
    assert s["length_delta"] == plus_sites[0]["length_delta"] == 1
    assert s["alt_protein"] == plus_sites[0]["alt_protein"]
    assert s["frameshift"] == plus_sites[0]["frameshift"]


# ---------------------------------------------------------------- 等价源 ID

def test_equivalent_source_ids_share_one_result():
    fasta = ">chr1\n" + "A" * 6 + CDS_REF + "\n"
    # 同聚物中的同一插入以两种写法出现 → 归并为一个标准位点
    rep = _report(fasta, csv_of(("a", 3, "A", "AA"), ("b", 5, "A", "AA")))
    assert len(rep["groups"]) == 1
    res = analyze_coding(rep["ref_seq"], rep["groups"],
                         [{"start": 7, "end": 18}], "+")
    assert len(res["sites"]) == 1
    assert res["sites"][0]["ids"] == ["a", "b"]


# ---------------------------------------------------------------- HTTP 流程

FASTA_APP = ">chr1\n" + CDS_REF + "\n"
CSV_APP = "ID,POS,REF,ALT\ns1,4,G,A\n"


@pytest.fixture()
def app(tmp_path):
    return create_app({
        "TESTING": True,
        "DATABASE": str(tmp_path / "t.db"),
        "WORKER_ENABLED": False,
        "COMPUTE_DELAY_SECONDS": 0,
    })


@pytest.fixture()
def client(app):
    return app.test_client()


def _upload(client, fasta=FASTA_APP, csv_text=CSV_APP):
    return client.post("/batches", data={
        "name": "t",
        "fasta": (io.BytesIO(fasta.encode()), "ref.fa"),
        "csv": (io.BytesIO(csv_text.encode()), "vars.csv"),
    }, follow_redirects=False)


def _run_coding(client, rid=1, start=1, end=12, strand="+", second=None):
    data = {"strand": strand, "exon1_start": str(start), "exon1_end": str(end)}
    if second:
        data["exon2_start"], data["exon2_end"] = map(str, second)
    return client.post(f"/batches/1/reports/{rid}/coding", data=data)


def test_http_analysis_view_shows_proteins(app, client):
    _upload(client)
    run_pending_jobs(app.config["DATABASE"])
    resp = _run_coding(client)
    assert resp.status_code == 302
    page = client.get(resp.headers["Location"]).data.decode()
    assert "MEF*" in page and "MKF*" in page
    assert "首个变化位置" in page and "编码可计算" in page
    assert "移码" in page and "提前终止" in page


def test_http_rejected_exon_definition_records_no_partial(app, client):
    _upload(client)
    run_pending_jobs(app.config["DATABASE"])
    resp = _run_coding(client, start=1, end=11)
    assert resp.status_code == 302
    conn = connect(app.config["DATABASE"])
    rows = conn.execute("SELECT * FROM analyses").fetchall()
    conn.close()
    assert len(rows) == 1
    assert rows[0]["status"] == "failed" and rows[0]["payload"] is None
    assert "3 的倍数" in rows[0]["error"]


def test_http_analysis_fixed_to_old_report_is_expired(app, client):
    _upload(client)
    run_pending_jobs(app.config["DATABASE"])
    _run_coding(client)                            # 报告 #1 的分析 #1
    client.post("/batches/1/inputs/csv", data={
        "file": (io.BytesIO(b"ID,POS,REF,ALT\nd1,4,GG,G\n"), "new.csv"),
    })
    run_pending_jobs(app.config["DATABASE"])

    page = client.get("/batches/1/reports/1/analyses/1").data.decode()
    assert "已过期" in page and "不是当前结论" in page
    assert "MEF*" in page                          # 旧依据仍可查
    report_page = client.get("/batches/1/reports/1").data.decode()
    assert "已过期" in report_page


def test_http_second_analysis_supersedes_first_within_same_report(app, client):
    _upload(client)
    run_pending_jobs(app.config["DATABASE"])
    _run_coding(client)
    _run_coding(client, start=1, end=12)
    first = client.get("/batches/1/reports/1/analyses/1").data.decode()
    second = client.get("/batches/1/reports/1/analyses/2").data.decode()
    assert "历史依据" in first and "当前结论" in second
    listing = client.get("/batches/1/reports/1").data.decode()
    assert "历史依据" in listing and "当前结论" in listing


def test_http_rejects_analysis_on_failed_report(app, client):
    _upload(client, csv_text="ID,POS,REF,ALT\nbad,4,C,G\n")
    run_pending_jobs(app.config["DATABASE"])
    resp = _run_coding(client)
    assert resp.status_code == 302
    conn = connect(app.config["DATABASE"])
    assert conn.execute("SELECT COUNT(*) c FROM analyses").fetchone()["c"] == 0
    conn.close()


def test_analysis_uses_report_fixed_reference_not_latest_input(app, client):
    """报告固定旧 FASTA；之后更正 FASTA，旧报告的分析仍按旧参考计算。"""
    _upload(client)
    run_pending_jobs(app.config["DATABASE"])
    # 新 FASTA 长度不变但序列不同（第 4 位改为 T）
    client.post("/batches/1/inputs/fasta", data={
        "file": (io.BytesIO((">chr1\n" + "ATGTAATTTTAA\n").encode()), "new.fa"),
    })
    run_pending_jobs(app.config["DATABASE"])

    resp = _run_coding(client, rid=1)               # 针对旧报告 #1
    page = client.get(resp.headers["Location"]).data.decode()
    assert "MEF*" in page and "MKF*" in page        # 依据旧参考 ATGGAATTTTAA
    assert "已过期" in page


def test_analysis_payload_not_leaked_into_failed(app, client):
    _upload(client)
    run_pending_jobs(app.config["DATABASE"])
    _run_coding(client, start=1, end=11)
    page = client.get("/batches/1/reports/1").data.decode()
    assert "已拒绝" in page
    detail = client.get("/batches/1/reports/1/analyses/1").data.decode()
    assert "整次拒绝原因" in detail and "3 的倍数" in detail
