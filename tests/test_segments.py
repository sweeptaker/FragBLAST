import math
import random

import pytest

from fragblast.segments import (
    columns_from_rows,
    find_fragments,
)


def _rows_to_columns(q_row: str, s_row: str):
    return columns_from_rows(q_row, s_row)


def _counts(columns):
    q_covered = [0]
    s_covered = [0]
    ident = [0]
    q = s = ids = 0
    for col in columns:
        if col.q_consumed:
            q += 1
        if col.s_consumed:
            s += 1
        if col.identical:
            ids += 1
        q_covered.append(q)
        s_covered.append(s)
        ident.append(ids)
    return q_covered, s_covered, ident


def brute_fragments(columns, query_len, subject_len, identity_pct, qcov_pct, scov_pct):
    """Reference implementation returning maximal ranges in residue space."""
    columns = tuple(columns)
    n = len(columns)
    q_covered, s_covered, ident = _counts(columns)
    q_pos = []
    s_pos = []
    q = s = -1
    for col in columns:
        q += 1 if col.q_consumed else 0
        s += 1 if col.s_consumed else 0
        q_pos.append(q)
        s_pos.append(s)
    min_q = math.floor(query_len * qcov_pct / 100) + 1
    min_s = math.floor(subject_len * scov_pct / 100) + 1

    valid = []
    for left in range(n):
        if not (columns[left].q_consumed and columns[left].s_consumed):
            continue
        for right in range(left, n):
            if not (columns[right].q_consumed and columns[right].s_consumed):
                continue
            length = right - left + 1
            ids = ident[right + 1] - ident[left]
            if ids * 100 <= identity_pct * length:
                continue
            q_len = q_covered[right + 1] - q_covered[left]
            s_len = s_covered[right + 1] - s_covered[left]
            if q_len * 100 <= qcov_pct * query_len:
                continue
            if s_len * 100 <= scov_pct * subject_len:
                continue
            valid.append((left, right))

    valid.sort(key=lambda lr: (lr[0], -lr[1]))
    maximal = []
    best_end = -1
    for left, right in valid:
        if right > best_end:
            best_end = right
            maximal.append((q_pos[left], q_pos[right]))
    return set(maximal)


def test_core_embedded_in_longer_low_identity_alignment():
    # Query/subject of length 100.  Alignment path covers 60 residues with 60%
    # identity overall (36/60), but the first 45 columns form a high-identity
    # core.  The tail is arranged all-mismatch so no trailing segment can
    # recover above threshold; brute force defines the expected maximal range.
    # First 45 columns: 36 identities (80%); next 15 columns: no identities.
    # Full 60-column alignment is 36/60 = 60% identity and therefore does not
    # qualify.  Brute force (against the same definition) gives the expected
    # maximal extension: identity stays above 70% while 6 of the 15 tail
    # columns are included, so the maximal fragment covers query residues 1-51.
    q_row = "A" * 36 + "B" * 24
    s_row = "A" * 36 + "C" * 24
    columns = _rows_to_columns(q_row, s_row)

    got = find_fragments(columns, 100, 100)
    expected = brute_fragments(columns, 100, 100, 70, 40, 40)
    got_ranges = {(f.q_start - 1, f.q_end - 1) for f in got}
    assert got_ranges == expected
    assert len(got) == 1
    fragment = got[0]
    assert (fragment.q_start, fragment.q_end) == (1, 51)
    assert fragment.identity > 70.0
    assert fragment.qcov > 40.0
    assert fragment.scov > 40.0
    # The 60-column full alignment itself must not qualify.
    assert fragment.aligned_len < 60


def test_two_islands_with_low_identity_linker_produce_staggered_maxima():
    # Two 41-residue islands (29/41 identical each) are separated by a
    # 10-residue all-mismatch linker.  Cross-linker ranges that include the
    # full pair would drop below the threshold, but the maximal ranges still
    # stagger: they may include the entire first island and the high-identity
    # head of the second island while staying above 70%.
    q_row = "A" * 29 + "B" * 12 + "B" * 10 + "A" * 29 + "B" * 12
    s_row = "A" * 29 + "C" * 12 + "D" * 10 + "A" * 29 + "C" * 12
    columns = _rows_to_columns(q_row, s_row)
    got = find_fragments(columns, 100, 100)
    got_ranges = {(f.q_start - 1, f.q_end - 1) for f in got}
    assert got_ranges == brute_fragments(columns, 100, 100, 70, 40, 40)
    assert len(got) > 2  # inclusion-maximal ranges may overlap/stagger
    for fragment in got:
        assert fragment.qcov > 40.0
        assert fragment.scov > 40.0
        assert fragment.identity > 70.0
    # Nothing is contained inside another reported range.
    by_start = sorted(got_ranges)
    kept = 0
    best_end = -1
    for start, end in by_start:
        if end > best_end:
            best_end = end
            kept += 1
    assert kept == len(got_ranges)


def test_internal_gap_counts_in_identity_but_not_coverage():
    # 45 core residue pairs + a 3-column internal query gap.  All 48 columns
    # contain 40 identities, so identity is 40/48, qcov is 45/100 and scov is
    # 48/100 (subject length 100).
    left = "A" * 22 + "B" * 2  # 22 identities over 24 columns
    gap = "C" * 3
    right = "A" * 18 + "B" * 3  # 18 identities over 21 columns
    q_row = left + "-" * 3 + right
    s_row = left + gap + right
    # "C" appears in both sides only within the gap columns, which for a query
    # gap have no query residue, so they are not identical residue pairs.
    s_row = "A" * 22 + "C" * 2 + "C" * 3 + "A" * 18 + "D" * 3
    q_row = "A" * 22 + "D" * 2 + "-" * 3 + "A" * 18 + "E" * 3
    columns = _rows_to_columns(q_row, s_row)
    fragments = find_fragments(columns, 100, 100, identity_pct=70, qcov_pct=40, scov_pct=40)
    assert len(fragments) == 1
    fragment = fragments[0]
    assert fragment.aligned_len == 48
    assert fragment.identities == 40
    assert math.isclose(fragment.identity, 100 * 40 / 48)
    assert math.isclose(fragment.qcov, 45.0)
    assert math.isclose(fragment.scov, 48.0)


def test_strict_threshold_at_100_returns_nothing():
    columns = _rows_to_columns("A" * 30, "A" * 30)
    assert find_fragments(columns, 30, 30, identity_pct=100, qcov_pct=0, scov_pct=0) == []
    assert len(find_fragments(columns, 30, 30, identity_pct=0, qcov_pct=0, scov_pct=0)) == 1


def test_random_paths_match_brute_force():
    rng = random.Random(20260907)
    letters = "ACDEFGHIKLMNPQRSTVWY"
    for case in range(80):
        n = rng.randint(8, 32)
        q_row_chars = []
        s_row_chars = []
        for _ in range(n):
            kind = rng.random()
            if kind < 0.08:
                q_row_chars.append("-")
                s_row_chars.append(rng.choice(letters))
            elif kind < 0.16:
                q_row_chars.append(rng.choice(letters))
                s_row_chars.append("-")
            else:
                qa = rng.choice(letters)
                sb = qa if rng.random() < 0.8 else rng.choice(letters)
                q_row_chars.append(qa)
                s_row_chars.append(sb)
        columns = columns_from_rows("".join(q_row_chars), "".join(s_row_chars))
        query_len = rng.randint(20, 70)
        subject_len = rng.randint(20, 70)
        ident_pct = rng.choice([50, 60, 70, 75, 85])
        cov_pct = rng.choice([30, 40, 45])
        expected = brute_fragments(
            columns, query_len, subject_len, ident_pct, cov_pct, cov_pct
        )
        got = find_fragments(
            columns,
            query_len,
            subject_len,
            identity_pct=ident_pct,
            qcov_pct=cov_pct,
            scov_pct=cov_pct,
        )
        got_ranges = {(f.q_start - 1, f.q_end - 1) for f in got}
        assert got_ranges == expected, f"case {case}: {got_ranges} != {expected}"


def test_rejects_bad_columns():
    with pytest.raises(ValueError):
        columns_from_rows("-", "-")
    with pytest.raises(ValueError):
        columns_from_rows("AC", "A")
