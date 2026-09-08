"""Thin wrapper around parasail's local (Smith-Waterman) protein alignment."""

from __future__ import annotations

from dataclasses import dataclass

from .segments import AlignmentColumn, columns_from_rows


@dataclass(frozen=True)
class PairAlignment:
    """A local alignment as a column vector plus 0-based residue offsets."""

    columns: tuple[AlignmentColumn, ...]
    q_offset: int
    s_offset: int
    score: int


def _import_parasail():
    try:
        import parasail
    except ImportError as exc:  # pragma: no cover - environment-specific
        raise RuntimeError(
            "parasail is required. Install it with: "
            "pip install parasail"
        ) from exc
    return parasail


def _matrix(parasail, name: str):
    try:
        return getattr(parasail, name)
    except AttributeError as exc:  # pragma: no cover - defensive
        raise ValueError(f"unknown scoring matrix: {name}") from exc


def _alignment_offsets(result, q_count: int, s_count: int) -> tuple[int, int]:
    """0-based positions of the first residues covered by the traceback.

    parasail reports the inclusive 0-based ``end_query``/``end_ref`` of the
    local alignment.  Walk back from those positions using the number of
    residues actually present in the traceback.  This is exact for local
    optimal alignments, which never begin or end with a gap column (trimming
    such a column would strictly improve the score).  The CIGAR object's
    ``beg_query``/``beg_ref`` fields are not reliable for this purpose.
    """
    q_end = getattr(result, "end_query", None)
    s_end = getattr(result, "end_ref", None)
    if q_end is None or s_end is None:
        raise RuntimeError(
            "cannot determine alignment start from parasail result"
        )
    return q_end - q_count + 1, s_end - s_count + 1


def align_pair(
    query: str,
    subject: str,
    *,
    gap_open: int = 11,
    gap_extend: int = 1,
    matrix: str = "blosum62",
) -> PairAlignment | None:
    """Return the optimal local alignment, or ``None`` when score <= 0."""
    if not query or not subject:
        raise ValueError("query and subject must be non-empty")
    parasail = _import_parasail()
    matrix_obj = _matrix(parasail, matrix)
    result = parasail.sw_trace_striped_sat(
        query, subject, gap_open, gap_extend, matrix_obj
    )
    traceback = result.traceback
    if result.score <= 0 or traceback is None:
        return None

    q_row = getattr(traceback, "query", None)
    s_row = getattr(traceback, "ref", None)
    if s_row is None:
        s_row = getattr(traceback, "subject", None)
    if q_row is None or s_row is None:
        raise RuntimeError(
            "unexpected parasail traceback object; cannot locate alignment rows"
        )

    columns = columns_from_rows(q_row, s_row)
    q_count = sum(1 for col in columns if col.q_consumed)
    s_count = sum(1 for col in columns if col.s_consumed)
    q_offset, s_offset = _alignment_offsets(result, q_count, s_count)
    return PairAlignment(
        columns=columns,
        q_offset=q_offset,
        s_offset=s_offset,
        score=int(result.score),
    )
