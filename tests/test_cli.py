from pathlib import Path

from fragblast.cli import main, read_fasta


def _write(path: Path, records):
    with path.open("w", encoding="utf-8") as handle:
        for ident, seq in records:
            handle.write(f">{ident} test record\n{seq}\n")


def test_read_fasta_handles_wrapped_and_blank_lines(tmp_path):
    fasta = tmp_path / "in.fa"
    _write(
        fasta,
        [
            ("q1", "MKTAYIAKQRQISFVKSHF"),
            ("s2", "ACD"),
        ],
    )
    assert read_fasta(str(fasta)) == [
        ("q1", "MKTAYIAKQRQISFVKSHF"),
        ("s2", "ACD"),
    ]


def test_cli_outputs_header_and_identical_hit(tmp_path, capsys):
    query = tmp_path / "q.fa"
    db = tmp_path / "db.fa"
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALP"
    _write(query, [("q1", seq)])
    _write(db, [("s1", seq)])
    out = tmp_path / "hits.tsv"
    rc = main(
        [
            "--query",
            str(query),
            "--db",
            str(db),
            "--out",
            str(out),
        ]
    )
    assert rc == 0
    lines = out.read_text(encoding="utf-8").strip().splitlines()
    assert lines[0].split("\t") == [
        "query_id",
        "subject_id",
        "q_start",
        "q_end",
        "s_start",
        "s_end",
        "aligned_len",
        "identity",
        "qcov",
        "scov",
    ]
    assert len(lines) == 2
    row = lines[1].split("\t")
    assert row[:6] == ["q1", "s1", "1", str(len(seq)), "1", str(len(seq))]
    assert float(row[7]) == 100.0


def _run_cli(tmp_path, query_seq, subject_seq, coverage_sequence):
    query = tmp_path / f"q_{coverage_sequence}.fa"
    db = tmp_path / f"db_{coverage_sequence}.fa"
    out = tmp_path / f"hits_{coverage_sequence}.tsv"
    _write(query, [("q1", query_seq)])
    _write(db, [("s1", subject_seq)])
    rc = main(
        [
            "--query",
            str(query),
            "--db",
            str(db),
            "--cov",
            "40",
            "--coverage_sequence",
            coverage_sequence,
            "--out",
            str(out),
        ]
    )
    assert rc == 0
    return out.read_text(encoding="utf-8").strip().splitlines()


def test_coverage_sequence_modes_select_which_coverage_is_enforced(tmp_path):
    core = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALP"
    long_tail = "W" * 200

    # Query is fully covered, but the subject hit covers only ~24% of subject.
    both_lines = _run_cli(tmp_path, core, core + long_tail, "both")
    query_lines = _run_cli(tmp_path, core, core + long_tail, "query")
    assert len(both_lines) == 1  # header only
    assert len(query_lines) == 2
    query_row = query_lines[1].split("\t")
    assert float(query_row[8]) > 40.0  # qcov
    assert float(query_row[9]) < 40.0  # scov

    # Same alignment, but now query carries the long tail and subject is fully
    # covered: only --coverage_sequence subject accepts the fragment.
    both_lines = _run_cli(tmp_path, core + long_tail, core, "both")
    subject_lines = _run_cli(tmp_path, core + long_tail, core, "subject")
    assert len(both_lines) == 1
    assert len(subject_lines) == 2
    subject_row = subject_lines[1].split("\t")
    assert float(subject_row[8]) < 40.0  # qcov
    assert float(subject_row[9]) > 40.0  # scov


def test_cli_rejects_empty_database(tmp_path):
    query = tmp_path / "q.fa"
    db = tmp_path / "db.fa"
    _write(query, [("q1", "MKTAYIAKQRQ")])
    _write(db, [])
    try:
        main(["--query", str(query), "--db", str(db)])
    except SystemExit as exc:
        assert exc.code != 0
    else:
        raise AssertionError("expected SystemExit for empty database")
