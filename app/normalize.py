"""变异解析、校验、左对齐标准化与等价归并。

纯函数模块：不依赖 Flask / SQLite，便于单独测试。
所有非法输入抛出 ValidationError，消息面向用户展示；调用方保证整批失败。
"""
from __future__ import annotations

import csv
import io

BASES = set("ACGTN")
MAX_REF_LEN = 4096
MAX_ROWS = 20

SNV = "SNV"
INS = "INS"
DEL = "DEL"

KIND_LABELS = {SNV: "单碱基替换", INS: "插入", DEL: "缺失"}


class ValidationError(Exception):
    """输入非法（文件级或行级）；消息面向用户展示。"""


# ---------------------------------------------------------------- 输入解析

def parse_fasta(text: str):
    """解析单条染色体 FASTA，返回 (chrom, seq)。序列统一大写。"""
    lines = [ln.strip() for ln in text.replace("\r\n", "\n").replace("\r", "\n").split("\n")]
    lines = [ln for ln in lines if ln]
    if not lines:
        raise ValidationError("FASTA 文件为空")
    if not lines[0].startswith(">"):
        raise ValidationError("FASTA 首行必须是以 '>' 开头的标题行")
    if any(ln.startswith(">") for ln in lines[1:]):
        raise ValidationError("只接受单条染色体 FASTA（发现多个标题行）")
    tokens = lines[0][1:].strip().split()
    if not tokens:
        raise ValidationError("FASTA 标题行缺少染色体名称")
    chrom = tokens[0]
    seq = "".join(lines[1:]).replace(" ", "").upper()
    if not seq:
        raise ValidationError("FASTA 序列为空")
    if len(seq) > MAX_REF_LEN:
        raise ValidationError(f"参考序列长度 {len(seq)} 超过上限 {MAX_REF_LEN}")
    bad = sorted(set(seq) - BASES)
    if bad:
        raise ValidationError(f"参考序列含非法字符 {''.join(bad)}（仅允许 A/C/G/T/N）")
    return chrom, seq


REQUIRED_COLS = ("id", "pos", "ref", "alt")


def parse_csv(text: str):
    """解析变异 CSV，返回行字典列表（id/pos/ref/alt/line）。"""
    try:
        rows = [r for r in csv.reader(io.StringIO(text))]
    except csv.Error as exc:
        raise ValidationError(f"CSV 解析失败: {exc}") from exc
    rows = [r for r in rows if any(cell.strip() for cell in r)]
    if not rows:
        raise ValidationError("CSV 文件为空")
    header = [c.strip().lower() for c in rows[0]]
    missing = [c for c in REQUIRED_COLS if c not in header]
    if missing:
        raise ValidationError("CSV 缺少必需列: " + ", ".join(c.upper() for c in missing))
    idx = {c: header.index(c) for c in REQUIRED_COLS}
    data = rows[1:]
    if not data:
        raise ValidationError("CSV 没有数据行")
    if len(data) > MAX_ROWS:
        raise ValidationError(f"CSV 数据行数 {len(data)} 超过上限 {MAX_ROWS}")

    out, seen = [], set()
    for n, r in enumerate(data, start=2):  # n = 含表头的文件行号
        def col(name, r=r):
            i = idx[name]
            return r[i].strip() if i < len(r) else ""

        vid = col("id")
        if not vid:
            raise ValidationError(f"第 {n} 行: ID 为空")
        if any(ch.isspace() for ch in vid) or "," in vid or ";" in vid:
            raise ValidationError(f"第 {n} 行: ID '{vid}' 含空白/逗号/分号")
        if vid in seen:
            raise ValidationError(f"第 {n} 行: ID '{vid}' 重复")
        seen.add(vid)

        pos_s = col("pos")
        if not pos_s.isdigit() or int(pos_s) < 1:
            raise ValidationError(f"第 {n} 行: POS '{pos_s}' 不是从 1 开始的正整数")
        pos = int(pos_s)

        ref, alt = col("ref").upper(), col("alt").upper()
        if not ref or not alt:
            raise ValidationError(f"第 {n} 行: REF/ALT 不能为空")
        for label, v in (("REF", ref), ("ALT", alt)):
            bad = sorted(set(v) - BASES)
            if bad:
                raise ValidationError(f"第 {n} 行: {label} 含非法字符 {''.join(bad)}")
        out.append({"id": vid, "pos": pos, "ref": ref, "alt": alt, "line": n})
    return out


# ---------------------------------------------------------------- 行级校验

def classify(row) -> str:
    """判定变异类型；不支持的写法抛 ValidationError。"""
    ref, alt, n = row["ref"], row["alt"], row["line"]
    if ref == alt:
        raise ValidationError(f"第 {n} 行: REF 与 ALT 相同")
    if len(ref) == 1 and len(alt) == 1:
        return SNV
    if len(ref) == 1 and len(alt) > 1 and alt.startswith(ref):
        return INS
    if len(alt) == 1 and len(ref) > 1 and ref.startswith(alt):
        return DEL
    raise ValidationError(
        f"第 {n} 行: 仅支持单碱基替换，或带锚定碱基的纯插入/缺失"
        f"（REF={ref}, ALT={alt}）"
    )


def check_ref(seq: str, row) -> None:
    """逐行核对 REF 与参考序列是否一致。"""
    pos, ref, n = row["pos"], row["ref"], row["line"]
    if pos > len(seq):
        raise ValidationError(f"第 {n} 行: POS {pos} 超出参考序列长度 {len(seq)}")
    end = pos - 1 + len(ref)
    if end > len(seq):
        raise ValidationError(
            f"第 {n} 行: REF 超出参考序列末端（POS {pos} 起 {len(ref)} 碱基，参考仅 {len(seq)}）"
        )
    actual = seq[pos - 1:end]
    if actual != ref:
        raise ValidationError(
            f"第 {n} 行: REF '{ref}' 与参考序列不符（该处参考为 '{actual}'）"
        )


# ---------------------------------------------------------------- 标准化

def normalize_variant(seq: str, row):
    """左对齐标准化，返回 (pos, ref, alt)（1 基 POS，含锚定碱基的 VCF 写法）。"""
    kind, pos, ref, alt = row["kind"], row["pos"], row["ref"], row["alt"]
    if kind == SNV:
        return pos, ref, alt
    if kind == INS:
        ins = alt[1:]
        p = pos - 1  # 锚定碱基的 0 基下标
        # 插入序列末碱基等于锚定碱基时，整体左移一位并旋转插入序列
        # （等价推导：seq[p] + S == S' + seq[p]  ⟺  S' = seq[p] + S[:-1]）
        while p >= 1 and seq[p] == ins[-1]:
            ins = seq[p] + ins[:-1]
            p -= 1
        return p + 1, seq[p], seq[p] + ins
    # DEL：s 为被删片段的 0 基起始下标（锚定在其左侧）
    dele = ref[1:]
    s = pos
    # 左移后锚定下标为 s-1，必须 >= 1（即 POS>=1 处仍有锚定碱基）
    while s >= 2 and seq[s - 1] == dele[-1]:
        dele = seq[s - 1] + dele[:-1]
        s -= 1
    return s, seq[s - 1] + dele, seq[s - 1]


def apply_variant(seq: str, row) -> str:
    """单条变异作用后的完整序列（等价判定依据）。"""
    i = row["pos"] - 1
    return seq[:i] + row["alt"] + seq[i + len(row["ref"]):]


# ---------------------------------------------------------------- 报告计算

def render_vcf(chrom: str, ref_len: int, groups) -> str:
    lines = [
        "##fileformat=VCFv4.2",
        "##source=variant-normalizer",
        f"##contig=<ID={chrom},length={ref_len}>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    for g in groups:
        lines.append("\t".join([
            chrom, str(g["pos"]), ";".join(g["ids"]), g["ref"], g["alt"], ".", ".", ".",
        ]))
    return "\n".join(lines) + "\n"


def compute_report(fasta_text: str, csv_text: str) -> dict:
    """整批计算：任一非法行 → ValidationError（整批失败，无部分结果）。"""
    chrom, seq = parse_fasta(fasta_text)
    rows = parse_csv(csv_text)
    for row in rows:  # 先全量校验，再计算，保证不产生半份结果
        row["kind"] = classify(row)
        check_ref(seq, row)

    groups = {}  # 完整突变序列 -> 归并组
    for row in rows:
        npos, nref, nalt = normalize_variant(seq, row)
        row["norm"] = (npos, nref, nalt)
        key = apply_variant(seq, row)  # 等价 = 单条变异作用后的完整序列相同
        g = groups.get(key)
        if g is None:
            g = {"pos": npos, "ref": nref, "alt": nalt, "kind": row["kind"], "ids": []}
            groups[key] = g
        g["ids"].append(row["id"])
        row["_key"] = key

    ordered = sorted(groups.values(), key=lambda g: (g["pos"], g["ref"], g["alt"]))
    for g in ordered:
        g["ids"].sort()
        g["kind_label"] = KIND_LABELS[g["kind"]]
    group_index = {id(g): i for i, g in enumerate(ordered)}

    out_rows = [{
        "id": row["id"], "line": row["line"],
        "pos": row["pos"], "ref": row["ref"], "alt": row["alt"],
        "kind": row["kind"], "kind_label": KIND_LABELS[row["kind"]],
        "norm_pos": row["norm"][0], "norm_ref": row["norm"][1], "norm_alt": row["norm"][2],
        "group": group_index[id(groups[row["_key"]])],
    } for row in rows]

    out_groups = [{k: g[k] for k in ("pos", "ref", "alt", "kind", "kind_label", "ids")}
                  for g in ordered]
    return {
        "chrom": chrom,
        "ref_len": len(seq),
        "groups": out_groups,
        "rows": out_rows,
        "vcf": render_vcf(chrom, len(seq), ordered),
    }
