"""编码影响分析：外显子拼接、CDS 重建、翻译与变异效应判定。

纯函数模块：不依赖 Flask / SQLite。调用方提供**报告固定版本**的参考序列与
已归并的标准位点（normalize.compute_report 的 groups），以及用户填写的外显子
定义。外显子定义非法时抛出 ValidationError，调用方必须整次拒绝，不得保留
部分位点结果。
"""
from __future__ import annotations

from .normalize import ValidationError

# 标准遗传密码（DNA，U 写作 T）
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

STATUS_CODING = "coding"
STATUS_NON_CODING = "non_coding"
STATUS_BOUNDARY = "boundary"

STATUS_LABELS = {
    STATUS_CODING: "编码可计算",
    STATUS_NON_CODING: "非编码",
    STATUS_BOUNDARY: "边界不可判定",
}

MAX_EXONS = 2


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def translate(seq: str) -> str:
    """按标准遗传密码翻译；遇到终止密码子即停止（含 '*'）。

    含 N 等未知碱基的密码子翻译为 'X'。末尾不足一个密码子的碱基忽略
    （参考 CDS 已保证长度为 3 的倍数；变异 CDS 可能移码）。
    """
    protein = []
    for i in range(0, len(seq) - 2, 3):
        aa = CODON_TABLE.get(seq[i:i + 3], "X")
        protein.append(aa)
        if aa == "*":
            break
    return "".join(protein)


# ---------------------------------------------------------------- 外显子校验

def validate_exons(ref_len: int, exons, strand: str):
    """校验外显子定义并返回按基因组坐标升序排列的列表。

    非法定义（数量、链、越界、起止颠倒、重叠、拼接长度非 3 倍数）整次拒绝。
    """
    if strand not in ("+", "-"):
        raise ValidationError("链方向必须为 +（正链）或 -（负链）")
    if not isinstance(exons, (list, tuple)) or not 1 <= len(exons) <= MAX_EXONS:
        raise ValidationError("必须提供 1 至 2 段编码外显子坐标")

    clean = []
    for i, exon in enumerate(exons, start=1):
        try:
            start, end = int(exon["start"]), int(exon["end"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValidationError(f"外显子 {i} 坐标必须是整数") from exc
        if start < 1 or end < 1 or end > ref_len or start > end:
            if start > end:
                raise ValidationError(
                    f"外显子 {i} 起点 {start} 大于终点 {end}")
            raise ValidationError(
                f"外显子 {i} 坐标 {start}-{end} 越界（参考序列长度 {ref_len}）")
        clean.append({"start": start, "end": end})

    clean.sort(key=lambda e: e["start"])
    if len(clean) == 2 and clean[0]["end"] >= clean[1]["start"]:
        raise ValidationError("两段外显子坐标重叠")

    total = sum(e["end"] - e["start"] + 1 for e in clean)
    if total % 3 != 0:
        raise ValidationError(
            f"拼接后外显子总长 {total} bp 不是 3 的倍数，无法按密码子翻译")
    return clean


# ---------------------------------------------------------------- 转录本构建

def _transcript(ref_seq: str, exons, strand: str):
    """按链方向拼接外显子，返回 (转录本顺序外显子, CDS)。"""
    ordered = sorted(exons, key=lambda e: e["start"], reverse=(strand == "-"))
    tx_exons, parts, offset = [], [], 0
    for exon in ordered:
        piece = ref_seq[exon["start"] - 1:exon["end"]]
        if strand == "-":
            piece = reverse_complement(piece)
        tx_exons.append({
            "start": exon["start"], "end": exon["end"],
            "offset": offset, "length": len(piece),
        })
        parts.append(piece)
        offset += len(piece)
    return tx_exons, "".join(parts)


def _cds_index(exon, pos: int, strand: str) -> int:
    """参考坐标 pos（1 基）在该外显子 CDS 中的 0 基下标。"""
    if strand == "+":
        return exon["offset"] + (pos - exon["start"])
    return exon["offset"] + (exon["end"] - pos)


def _containing_exon(tx_exons, start: int, end: int):
    """返回完整覆盖 [start, end] 的外显子；没有则返回 None。"""
    for exon in tx_exons:
        if exon["start"] <= start and end <= exon["end"]:
            return exon
    return None


def _overlaps_any(tx_exons, start: int, end: int) -> bool:
    return any(start <= e["end"] and e["start"] <= end for e in tx_exons)


# ---------------------------------------------------------------- 位点重建

def _rebuild(ref_seq, tx_exons, cds: str, strand: str, group):
    """尝试重建某标准位点的变异 CDS。

    返回 (status, alt_cds, reason, detail)。只有编码可计算时 alt_cds 非空。
    规则：插入点与受影响碱基都必须完整落在**同一段**外显子内；跨外显子
    边界或紧贴边界 → 边界不可判定；完全不涉及外显子 → 非编码。
    """
    pos, ref, alt, kind = group["pos"], group["ref"], group["alt"], group["kind"]

    if kind == "SNV":
        exon = _containing_exon(tx_exons, pos, pos)
        if exon is None:
            return (STATUS_NON_CODING, None,
                    "标准位点不在任何编码外显子内", {})
        idx = _cds_index(exon, pos, strand)
        base = alt if strand == "+" else reverse_complement(alt)
        alt_cds = cds[:idx] + base + cds[idx + 1:]
        return (STATUS_CODING, alt_cds, None,
                {"cds_start": idx + 1, "cds_end": idx + 1})

    if kind == "INS":
        # 插入点位于 POS 与 POS+1 之间；两侧参考碱基都必须在外显子内部。
        if pos >= len(ref_seq):
            exon = _containing_exon(tx_exons, pos, pos)
            if exon is not None:
                return (STATUS_BOUNDARY, None,
                        "插入点位于参考序列末端，超出可判定的编码边界", {})
            return (STATUS_NON_CODING, None,
                    "插入点不在任何编码外显子内", {})
        exon = _containing_exon(tx_exons, pos, pos + 1)
        if exon is None:
            if _overlaps_any(tx_exons, pos, pos + 1):
                return (STATUS_BOUNDARY, None,
                        "插入点位于外显子边界或跨外显子，无法判定编码影响", {})
            return (STATUS_NON_CODING, None,
                    "插入点不在任何编码外显子内", {})
        cut = _cds_index(exon, pos, strand) + (1 if strand == "+" else 0)
        inserted = alt[1:]
        if strand == "-":
            inserted = reverse_complement(inserted)
        alt_cds = cds[:cut] + inserted + cds[cut:]
        return (STATUS_CODING, alt_cds, None,
                {"cds_start": cut + 1, "cds_end": cut + len(inserted)})

    # DEL：REF 带一个左锚定碱基，真正受影响的是 POS+1 .. POS+len(REF)-1
    start, end = pos + 1, pos + len(ref) - 1
    exon = _containing_exon(tx_exons, start, end)
    if exon is not None:
        length = end - start + 1
        if strand == "+":
            cds_start = _cds_index(exon, start, strand)
        else:
            cds_start = _cds_index(exon, end, strand)
        alt_cds = cds[:cds_start] + cds[cds_start + length:]
        return (STATUS_CODING, alt_cds, None,
                {"cds_start": cds_start + 1, "cds_end": cds_start + length})
    if _overlaps_any(tx_exons, start, end):
        return (STATUS_BOUNDARY, None,
                "缺失区间跨外显子边界或跨越两段外显子，无法判定编码影响", {})
    return (STATUS_NON_CODING, None,
            "受影响碱基不在任何编码外显子内", {})


def first_changed_position(ref_protein: str, alt_protein: str):
    """首个变化的氨基酸位置（1 基）；完全相同（或无氨基酸变化）返回 None。"""
    for i, (a, b) in enumerate(zip(ref_protein, alt_protein)):
        if a != b:
            return i + 1
    if len(ref_protein) != len(alt_protein):
        return min(len(ref_protein), len(alt_protein)) + 1
    return None


def _first_stop(protein: str):
    idx = protein.find("*")
    return idx if idx >= 0 else None


def is_premature_stop(ref_protein: str, alt_protein: str) -> bool:
    """变异蛋白的终止密码子早于参考蛋白的终止密码子。"""
    alt_stop = _first_stop(alt_protein)
    if alt_stop is None:
        return False
    ref_stop = _first_stop(ref_protein)
    return ref_stop is None or alt_stop < ref_stop


# ---------------------------------------------------------------- 入口

def analyze_coding(ref_seq: str, groups, exons, strand: str) -> dict:
    """对报告中每个标准位点重建变异 CDS 并翻译。

    外显子非法 → ValidationError（整次拒绝）。位点不在外显子内或处于边界
    不是失败：逐位点标记为非编码 / 边界不可判定。
    """
    exons = validate_exons(len(ref_seq), exons, strand)
    tx_exons, cds = _transcript(ref_seq, exons, strand)
    ref_protein = translate(cds)

    sites, counts = [], {STATUS_CODING: 0, STATUS_NON_CODING: 0, STATUS_BOUNDARY: 0}
    for i, group in enumerate(groups, start=1):
        status, alt_cds, reason, detail = _rebuild(
            ref_seq, tx_exons, cds, strand, group)
        counts[status] += 1
        site = {
            "index": i,
            "ids": list(group.get("ids", [])),
            "pos": group["pos"],
            "ref": group["ref"],
            "alt": group["alt"],
            "kind": group["kind"],
            "kind_label": group.get("kind_label", group["kind"]),
            "status": status,
            "status_label": STATUS_LABELS[status],
            "reason": reason,
            "frameshift": False,
            "premature_stop": False,
            "first_change": None,
            "ref_protein": None,
            "alt_protein": None,
            "length_delta": 0,
        }
        if status == STATUS_CODING:
            alt_protein = translate(alt_cds)
            delta = len(alt_cds) - len(cds)
            site.update({
                "ref_protein": ref_protein,
                "alt_protein": alt_protein,
                "first_change": first_changed_position(ref_protein, alt_protein),
                "frameshift": delta % 3 != 0,
                "premature_stop": is_premature_stop(ref_protein, alt_protein),
                "length_delta": delta,
                "cds_start": detail.get("cds_start"),
                "cds_end": detail.get("cds_end"),
            })
        sites.append(site)

    return {
        "exons": exons,
        "strand": strand,
        "transcript_exons": tx_exons,
        "cds_len": len(cds),
        "cds_seq": cds,
        "ref_protein": ref_protein,
        "sites": sites,
        "summary": counts,
    }
