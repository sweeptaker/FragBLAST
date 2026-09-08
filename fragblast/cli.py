"""Command-line entry point for FragBLAST."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .align import align_pair
from .segments import find_fragments


def read_fasta(path: str) -> list[tuple[str, str]]:
    """Parse a FASTA file; header id is the first whitespace-delimited token."""
    sequences: list[tuple[str, str]] = []
    current_id: str | None = None
    chunks: list[str] = []
    try:
        handle = open(path, encoding="utf-8")
    except OSError as exc:
        raise ValueError(f"cannot open {path!r}: {exc}") from exc
    with handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences.append((current_id, "".join(chunks)))
                fields = line[1:].split()
                if not fields:
                    raise ValueError(f"empty FASTA header at {path!r}")
                current_id = fields[0]
                chunks = []
            else:
                if current_id is None:
                    raise ValueError(
                        f"sequence data before first header in {path!r}"
                    )
                chunks.append(line)
    if current_id is not None:
        sequences.append((current_id, "".join(chunks)))
    return sequences


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m fragblast",
        description=(
            "Re-align each query against each subject (parasail local "
            "alignment), then report maximal fragments whose identity, query "
            "coverage and subject coverage exceed the given thresholds."
        ),
    )
    parser.add_argument("--query", required=True, help="query proteins (FASTA)")
    parser.add_argument("--db", required=True, help="subject proteins (FASTA)")
    parser.add_argument("--out", help="output TSV path (default: stdout)")
    parser.add_argument(
        "--ident",
        type=float,
        default=70.0,
        help="minimum percent identity, strict (default: 70)",
    )
    parser.add_argument(
        "--cov",
        type=float,
        default=40.0,
        help="minimum coverage percent, strict (default: 40)",
    )
    parser.add_argument(
        "--coverage_sequence",
        choices=("both", "query", "subject"),
        default="both",
        help=(
            "which sequences the --cov threshold applies to: both query and "
            "subject (default), query only, or subject only"
        ),
    )
    parser.add_argument(
        "--gap-open",
        type=int,
        default=11,
        help="gap open penalty (default: 11)",
    )
    parser.add_argument(
        "--gap-extend",
        type=int,
        default=1,
        help="gap extend penalty (default: 1)",
    )
    parser.add_argument(
        "--matrix",
        default="blosum62",
        help="parasail scoring matrix name (default: blosum62)",
    )
    return parser


def _write_row(writer, values) -> None:
    writer.write("\t".join(str(v) for v in values))
    writer.write("\n")


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        queries = read_fasta(args.query)
        subjects = read_fasta(args.db)
    except ValueError as exc:
        parser.error(str(exc))

    if not queries:
        parser.error(f"no sequences found in query file {args.query!r}")
    if not subjects:
        parser.error(f"no sequences found in db file {args.db!r}")
    if not 0 <= args.ident <= 100:
        parser.error("--ident must be between 0 and 100")
    if not 0 <= args.cov <= 100:
        parser.error("--cov must be between 0 and 100")

    if args.coverage_sequence in ("both", "query"):
        qcov_pct = args.cov
    else:
        qcov_pct = 0.0
    if args.coverage_sequence in ("both", "subject"):
        scov_pct = args.cov
    else:
        scov_pct = 0.0

    handle = open(args.out, "w", encoding="utf-8") if args.out else sys.stdout
    header = [
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
    with handle as writer:
        _write_row(writer, header)
        for query_id, query_seq in queries:
            for subject_id, subject_seq in subjects:
                alignment = align_pair(
                    query_seq,
                    subject_seq,
                    gap_open=args.gap_open,
                    gap_extend=args.gap_extend,
                    matrix=args.matrix,
                )
                if alignment is None:
                    continue
                fragments = find_fragments(
                    alignment.columns,
                    len(query_seq),
                    len(subject_seq),
                    identity_pct=args.ident,
                    qcov_pct=qcov_pct,
                    scov_pct=scov_pct,
                    q_offset=alignment.q_offset,
                    s_offset=alignment.s_offset,
                )
                for frag in fragments:
                    _write_row(
                        writer,
                        (
                            query_id,
                            subject_id,
                            frag.q_start,
                            frag.q_end,
                            frag.s_start,
                            frag.s_end,
                            frag.aligned_len,
                            f"{frag.identity:.6g}",
                            f"{frag.qcov:.6g}",
                            f"{frag.scov:.6g}",
                        ),
                    )
    return 0
