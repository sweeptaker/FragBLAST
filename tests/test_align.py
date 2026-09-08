import pytest

from fragblast.align import align_pair


def test_identical_sequence_aligns_over_full_length():
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALP"
    alignment = align_pair(seq, seq)
    assert alignment is not None
    assert alignment.score > 0
    assert alignment.q_offset == 0
    assert alignment.s_offset == 0
    assert len(alignment.columns) == len(seq)
    assert all(col.q_consumed and col.s_consumed and col.identical for col in alignment.columns)


def test_local_alignment_offsets_are_recovered():
    core = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALP"
    query = "A" * 15 + core + "W" * 15
    alignment = align_pair(query, core)
    assert alignment is not None
    assert alignment.q_offset == 15
    assert alignment.s_offset == 0
    assert len(alignment.columns) == len(core)
    assert all(col.identical for col in alignment.columns)


def test_divergent_sequences_may_produce_no_alignment():
    a = "A" * 30
    b = "W" * 30
    assert align_pair(a, b) is None


def test_unknown_matrix_is_rejected():
    with pytest.raises(ValueError):
        align_pair("ACDEFGHIK", "ACDEFGHIK", matrix="not-a-matrix")
