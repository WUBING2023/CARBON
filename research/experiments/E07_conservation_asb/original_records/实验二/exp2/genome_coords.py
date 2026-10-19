from __future__ import annotations


def norm_chrom(c: str) -> str:
    return str(c).replace("chr", "")


def array_index_to_genomic(cs: int, ce: int, strand: str, idx: int) -> int:
    st = str(strand).strip()
    if st == "+":
        return int(cs) + int(idx)
    if st == "-":
        return int(ce) - 1 - int(idx)
    raise ValueError(f"bad strand: {strand!r}")


def genomic_to_array_index(cs: int, ce: int, strand: str, g: int) -> int | None:
    if g < int(cs) or g >= int(ce):
        return None
    st = str(strand).strip()
    if st == "+":
        return int(g) - int(cs)
    if st == "-":
        return int(ce) - 1 - int(g)
    raise ValueError(f"bad strand: {strand!r}")


def seqlet_genomic_span(
    cs: int,
    ce: int,
    strand: str,
    start: int,
    end: int,
) -> tuple[int, int]:
    """seqlet 在窗内 [start, end) 对应基因组闭区间 [g_lo, g_hi]。"""
    g0 = array_index_to_genomic(cs, ce, strand, start)
    g1 = array_index_to_genomic(cs, ce, strand, end - 1)
    return (min(g0, g1), max(g0, g1))


def seqlet_overlaps_half_open_interval(
    g_lo: int,
    g_hi: int,
    interval_start: int,
    interval_end_excl: int,
) -> bool:
    """g_lo..g_hi 闭区间与 [interval_start, interval_end_excl) 是否有交集。"""
    return g_lo < int(interval_end_excl) and g_hi >= int(interval_start)
