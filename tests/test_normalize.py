"""归一化核心逻辑测试：左对齐、等价归并、逐行 REF 校验、整批失败。"""
import pytest

from app.normalize import (ValidationError, apply_variant, compute_report,
                           parse_csv, parse_fasta)

FASTA = ">chr1\n" + "A" * 20 + "\n"  # 纯 A 序列，便于构造同聚物场景


def csv_of(*rows):
    return "ID,POS,REF,ALT\n" + "\n".join(",".join(map(str, r)) for r in rows) + "\n"


# ---------------------------------------------------------------- 输入解析

def test_fasta_ok_multiline_and_lowercase():
    chrom, seq = parse_fasta(">chrX some description\nacgt\nACGT\n")
    assert chrom == "chrX"
    assert seq == "ACGTACGT"


def test_fasta_rejects_multiple_sequences():
    with pytest.raises(ValidationError, match="单条染色体"):
        parse_fasta(">a\nACGT\n>b\nACGT\n")


def test_fasta_rejects_over_4096():
    with pytest.raises(ValidationError, match="4096"):
        parse_fasta(">c\n" + "A" * 4097 + "\n")


def test_fasta_accepts_exactly_4096():
    _, seq = parse_fasta(">c\n" + "A" * 4096 + "\n")
    assert len(seq) == 4096


def test_fasta_rejects_bad_chars():
    with pytest.raises(ValidationError, match="非法字符"):
        parse_fasta(">c\nACGTX\n")


def test_csv_rejects_over_20_rows():
    rows = "\n".join(f"v{i},1,A,C" for i in range(21))
    with pytest.raises(ValidationError, match="20"):
        parse_csv(f"ID,POS,REF,ALT\n{rows}\n")


def test_csv_rejects_duplicate_ids():
    with pytest.raises(ValidationError, match="重复"):
        parse_csv(csv_of(("v1", 1, "A", "C"), ("v1", 2, "A", "C")))


def test_csv_rejects_pos_zero_and_non_numeric():
    with pytest.raises(ValidationError, match="正整数"):
        parse_csv(csv_of(("v1", 0, "A", "C")))
    with pytest.raises(ValidationError, match="正整数"):
        parse_csv(csv_of(("v1", "x", "A", "C")))


def test_csv_header_case_insensitive_and_reordered():
    rows = parse_csv("pos,ALT,ref,ID\n1,c,a,v1\n")
    assert rows[0]["ref"] == "A" and rows[0]["alt"] == "C" and rows[0]["id"] == "v1"


def test_csv_missing_column_fails():
    with pytest.raises(ValidationError, match="缺少必需列"):
        parse_csv("ID,POS,REF\nv1,1,A\n")


# ---------------------------------------------------------------- 行级校验

def test_ref_mismatch_fails_whole_batch():
    csv_text = csv_of(("good", 1, "A", "C"), ("bad", 5, "C", "G"))
    with pytest.raises(ValidationError, match="与参考序列不符"):
        compute_report(FASTA, csv_text)


def test_ref_equals_alt_fails():
    with pytest.raises(ValidationError, match="相同"):
        compute_report(FASTA, csv_of(("v1", 3, "A", "A")))


def test_mnp_not_supported():
    with pytest.raises(ValidationError, match="纯插入/缺失"):
        compute_report(FASTA, csv_of(("v1", 1, "AA", "CT")))


def test_indel_without_anchor_not_supported():
    # REF/ALT 首碱基不一致的"插入"不是带锚定碱基的纯插入
    with pytest.raises(ValidationError, match="纯插入/缺失"):
        compute_report(FASTA, csv_of(("v1", 1, "A", "CA")))


def test_pos_beyond_reference_fails():
    with pytest.raises(ValidationError, match="超出参考序列"):
        compute_report(FASTA, csv_of(("v1", 21, "A", "C")))


def test_deletion_overhang_fails():
    with pytest.raises(ValidationError, match="超出参考序列末端"):
        compute_report(FASTA, csv_of(("v1", 20, "AA", "A")))


# ---------------------------------------------------------------- 左对齐

def test_insertion_left_aligned_in_homopolymer():
    rep = compute_report(FASTA, csv_of(("ins", 10, "A", "AA")))
    g = rep["groups"][0]
    assert (g["pos"], g["ref"], g["alt"]) == (1, "A", "AA")


def test_deletion_left_aligned_in_homopolymer():
    rep = compute_report(FASTA, csv_of(("del", 10, "AA", "A")))
    g = rep["groups"][0]
    assert (g["pos"], g["ref"], g["alt"]) == (1, "AA", "A")


def test_deletion_left_aligned_in_repeat_with_rotation():
    fasta = ">c\nTGCAGCAGT\n"
    rep = compute_report(fasta, csv_of(("d", 2, "GCAG", "G")))  # 删 0 基 2..4 的 CAG
    g = rep["groups"][0]
    assert (g["pos"], g["ref"], g["alt"]) == (1, "TGCA", "T")


def test_insertion_rotation_in_tandem_repeat():
    fasta = ">c\nGATATG\n"  # G A T A T G
    # 在 pos3(T) 后插 AT 与在 pos4(A) 后插 TA 是同一变异（突变序列均为 GATATATG）
    rep = compute_report(fasta, csv_of(("a", 3, "T", "TAT"), ("b", 4, "A", "ATA")))
    assert len(rep["groups"]) == 1
    g = rep["groups"][0]
    # 最左可表示位置：POS=1，插入序列旋转为 AT
    assert (g["pos"], g["ref"], g["alt"]) == (1, "G", "GAT")
    assert g["ids"] == ["a", "b"]


def test_deletion_rotation_in_tandem_repeat():
    fasta = ">c\nGATATG\n"  # G A T A T G
    # 删 pos2-3 的 AT 与删 pos4-5 的 AT 是同一变异（突变序列均为 GATG）
    rep = compute_report(fasta, csv_of(("x", 1, "GAT", "G"), ("y", 3, "TAT", "T")))
    assert len(rep["groups"]) == 1
    g = rep["groups"][0]
    assert (g["pos"], g["ref"], g["alt"]) == (1, "GAT", "G")
    assert g["ids"] == ["x", "y"]


def test_leftmost_insertion_at_position_1():
    rep = compute_report(FASTA, csv_of(("ins", 1, "A", "AA")))
    assert rep["groups"][0]["pos"] == 1


# ---------------------------------------------------------------- 等价归并

def test_equivalent_indels_merged_with_all_ids():
    csv_text = csv_of(
        ("i1", 3, "A", "AA"),      # 同聚物中不同位置的同一插入
        ("i2", 7, "A", "AA"),
        ("d1", 12, "AA", "A"),     # 同一缺失的两种写法
        ("d2", 15, "AA", "A"),
        ("s1", 5, "A", "C"),       # 独立 SNV
    )
    rep = compute_report(FASTA, csv_text)
    assert len(rep["groups"]) == 3
    by_kind = {g["kind"]: g for g in rep["groups"]}
    assert by_kind["INS"]["ids"] == ["i1", "i2"]
    assert by_kind["INS"]["pos"] == 1
    assert by_kind["DEL"]["ids"] == ["d1", "d2"]
    assert by_kind["DEL"]["pos"] == 1
    assert by_kind["SNV"]["ids"] == ["s1"]
    assert (by_kind["SNV"]["pos"], by_kind["SNV"]["alt"]) == (5, "C")


def test_duplicate_snv_rows_merged():
    rep = compute_report(FASTA, csv_of(("a", 4, "A", "T"), ("b", 4, "A", "T")))
    assert len(rep["groups"]) == 1
    assert rep["groups"][0]["ids"] == ["a", "b"]


def test_distinct_variants_not_merged():
    rep = compute_report(FASTA, csv_of(("a", 4, "A", "T"), ("b", 4, "A", "C")))
    assert len(rep["groups"]) == 2


def test_apply_variant_uses_full_sequence_semantics():
    fasta = ">c\nACGT\n"
    rows = parse_csv(csv_of(("v", 2, "C", "T")))
    assert apply_variant("ACGT", rows[0]) == "ATGT"


# ---------------------------------------------------------------- VCF 输出

def test_vcf_output_standard_format():
    rep = compute_report(FASTA, csv_of(("i1", 3, "A", "AA"), ("i2", 9, "A", "AA")))
    lines = rep["vcf"].splitlines()
    assert lines[0] == "##fileformat=VCFv4.2"
    assert lines[2] == "##contig=<ID=chr1,length=20>"
    assert lines[3].startswith("#CHROM")
    assert lines[4] == "chr1\t1\ti1;i2\tA\tAA\t.\t.\t."


def test_rows_mapping_present():
    rep = compute_report(FASTA, csv_of(("i1", 3, "A", "AA")))
    row = rep["rows"][0]
    assert row["id"] == "i1" and row["pos"] == 3
    assert (row["norm_pos"], row["norm_ref"], row["norm_alt"]) == (1, "A", "AA")
    assert row["group"] == 0
