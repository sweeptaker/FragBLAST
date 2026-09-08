"""Maximal fragment enumeration along one gapped alignment path.

A fragment is a contiguous range of alignment columns.  It qualifies when

* percent identity (identical residue columns / all columns in the range,
  internal gaps count as non-identical) is strictly greater than
  ``identity_pct``;
* the number of distinct query residues covered is strictly greater than
  ``qcov_pct`` percent of the full query length; and
* the same holds for subject residues against the full subject length.

Only inclusion-maximal qualifying fragments are returned (no qualifying
sub-fragment contained in another qualifying fragment is reported).
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


@dataclass(frozen=True)
class AlignmentColumn:
    """One column of a pairwise alignment.

    ``q_consumed`` / ``s_consumed`` tell whether this column contains a
    residue from the query / subject (a gap is the opposite side).
    ``identical`` is true when both sides contain the same residue.
    """

    q_consumed: bool
    s_consumed: bool
    identical: bool


@dataclass(frozen=True)
class Fragment:
    """A maximal qualifying fragment, with 1-based inclusive coordinates."""

    q_start: int
    q_end: int
    s_start: int
    s_end: int
    aligned_len: int
    identities: int
    identity: float  # percent
    qcov: float  # percent of full query length
    scov: float  # percent of full subject length


def columns_from_rows(q_row: str, s_row: str) -> tuple[AlignmentColumn, ...]:
    """Turn two equal-length traceback rows (gaps shown as ``-``) into columns."""
    if len(q_row) != len(s_row):
        raise ValueError("alignment rows have different lengths")
    columns = []
    for q_ch, s_ch in zip(q_row, s_row):
        q_ok = q_ch != "-"
        s_ok = s_ch != "-"
        if not q_ok and not s_ok:
            raise ValueError("alignment contains an empty column")
        columns.append(
            AlignmentColumn(
                q_consumed=q_ok,
                s_consumed=s_ok,
                identical=q_ok and s_ok and q_ch == s_ch,
            )
        )
    return tuple(columns)


def _pct_fraction(value: float) -> Fraction:
    return Fraction(str(value)) / 100


def _rightmost_greater(values: list[int], ql: int, qr: int, x: int) -> int:
    """Return the largest index in ``[ql, qr]`` with ``values[i] > x``."""
    if ql > qr:
        return -1
    n = len(values)
    size = 1
    while size < n:
        size *= 2
    tree = [-(1 << 62)] * (2 * size)
    tree[size : size + n] = values
    for i in range(size - 1, 0, -1):
        tree[i] = max(tree[2 * i], tree[2 * i + 1])

    # Depth-first search, right child first, so the first positive leaf is the
    # rightmost one.  Nodes whose whole range cannot contain a qualifying value
    # are pruned.
    stack = [(1, 0, size - 1)]
    while stack:
        node, lo, hi = stack.pop()
        if hi < ql or lo > qr or tree[node] <= x:
            continue
        if lo == hi:
            return lo if lo < n else -1
        mid = (lo + hi) // 2
        # Push left first so that the right child is explored first.
        stack.append((node * 2, lo, mid))
        stack.append((node * 2 + 1, mid + 1, hi))
    return -1


def find_fragments(
    columns,
    query_len: int,
    subject_len: int,
    identity_pct: float = 70.0,
    qcov_pct: float = 40.0,
    scov_pct: float = 40.0,
    *,
    q_offset: int = 0,
    s_offset: int = 0,
) -> list[Fragment]:
    """Enumerate inclusion-maximal qualifying fragments on one alignment path."""
    cols = tuple(columns)
    if not cols:
        return []
    if query_len <= 0 or subject_len <= 0:
        raise ValueError("query_len and subject_len must be positive")
    if q_offset < 0 or s_offset < 0:
        raise ValueError("offsets must be non-negative")

    t_ident = _pct_fraction(identity_pct)
    t_q = _pct_fraction(qcov_pct)
    t_s = _pct_fraction(scov_pct)
    if any(t < 0 or t > 1 for t in (t_ident, t_q, t_s)):
        raise ValueError("threshold percentages must be within [0, 100]")

    # Strictly-greater thresholds: the smallest integer residue count that
    # exceeds ``length * t``.
    min_q = (query_len * t_q.numerator) // t_q.denominator + 1
    min_s = (subject_len * t_s.numerator) // t_s.denominator + 1
    if min_q > query_len or min_s > subject_len:
        return []

    n = len(cols)
    den = t_ident.denominator
    num = t_ident.numerator
    w_ident = den - num
    w_other = -num

    # Prefix arrays.  P is the identity-threshold weighted prefix sum:
    #   P[i] = sum_{k<i} (1-t_ident if column k is identical else -t_ident),
    # scaled by ``den`` so all arithmetic stays integral.  A range [l, r) is
    # valid w.r.t. identity iff P[r] - P[l] > 0.
    prefix = [0] * (n + 1)
    q_covered = [0] * (n + 1)
    s_covered = [0] * (n + 1)
    ident_count = [0] * (n + 1)
    q_pos = [-1] * n
    s_pos = [-1] * n

    q_so_far = 0
    s_so_far = 0
    for i, col in enumerate(cols):
        prefix[i + 1] = prefix[i] + (w_ident if col.identical else w_other)
        if col.q_consumed:
            q_pos[i] = q_so_far
            q_so_far += 1
        if col.s_consumed:
            s_pos[i] = s_so_far
            s_so_far += 1
        q_covered[i + 1] = q_so_far
        s_covered[i + 1] = s_so_far
        ident_count[i + 1] = ident_count[i] + (1 if col.identical else 0)

    # Canonical boundary columns contain a residue from both sequences.
    endpoints = [i for i in range(n) if q_pos[i] >= 0 and s_pos[i] >= 0]
    if not endpoints:
        return []

    endpoint_index = [-1] * n
    for order, col_idx in enumerate(endpoints):
        endpoint_index[col_idx] = order

    # Prefix sum value just after each canonical column (i.e. right boundary).
    endpoint_prefix = [prefix[col_idx + 1] for col_idx in endpoints]
    m = len(endpoints)

    # For every possible left boundary l, the only candidate that can be
    # inclusion-maximal is [l, r_last] where r_last is the *last* canonical
    # column whose cumulative weighted sum remains positive.  Coverage only
    # grows as the range grows, so checking r_last is sufficient.
    candidates = []
    for l_idx in endpoints:
        pos_l = endpoint_index[l_idx]
        best_order = _rightmost_greater(
            endpoint_prefix, pos_l, m - 1, prefix[l_idx]
        )
        if best_order < 0:
            continue
        r_idx = endpoints[best_order]
        q_span = q_covered[r_idx + 1] - q_covered[l_idx]
        s_span = s_covered[r_idx + 1] - s_covered[l_idx]
        if q_span >= min_q and s_span >= min_s:
            candidates.append((l_idx, r_idx))

    # Remove candidates contained in another candidate.  Sorting by start
    # ascending then end descending makes a previous range a superset exactly
    # when its end is at least as large.
    candidates.sort(key=lambda lr: (lr[0], -lr[1]))
    fragments = []
    best_end = -1
    for l_idx, r_idx in candidates:
        if r_idx <= best_end:
            continue
        best_end = r_idx

        aligned_len = r_idx - l_idx + 1
        ids = ident_count[r_idx + 1] - ident_count[l_idx]
        q_len = q_covered[r_idx + 1] - q_covered[l_idx]
        s_len = s_covered[r_idx + 1] - s_covered[l_idx]
        fragments.append(
            Fragment(
                q_start=q_offset + q_pos[l_idx] + 1,
                q_end=q_offset + q_pos[r_idx] + 1,
                s_start=s_offset + s_pos[l_idx] + 1,
                s_end=s_offset + s_pos[r_idx] + 1,
                aligned_len=aligned_len,
                identities=ids,
                identity=100.0 * ids / aligned_len,
                qcov=100.0 * q_len / query_len,
                scov=100.0 * s_len / subject_len,
            )
        )
    return fragments
