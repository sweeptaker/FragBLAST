#!/usr/bin/env python
"""DIAMOND + FragBLAST functional-fragment discovery pipeline.

1. DIAMOND finds candidate query-subject pairs.
2. Hard-threshold hits are reported as "solid" functional genes (ARGs in the
   accompanying case study).
3. Remaining pairs that satisfy two necessary conditions are treated as
   "potential":
   - coverage of the relevant sequence(s) is already above ``--cov``;
   - ``pident * coverage > ident * cov`` (a necessary condition for any
     embedded fragment that passes both hard thresholds).
4. FragBLAST re-aligns each potential query only against its matching
   subjects (the cheapest candidate reduction option) and enumerates maximal
   qualifying fragments along the local alignment path.
5. Queries with at least one qualifying fragment are reported as "new"
   functional genes.

All application thresholds are strictly greater-than. Coverage is controlled
by ``--coverage_sequence`` (query, subject, or both), and DIAMOND pre-filters
are derived automatically from the requested thresholds instead of being set
by hand.
"""

from __future__ import annotations

import argparse
import gzip
import math
import multiprocessing as mp
import shutil
import subprocess
import sys
from pathlib import Path

# Make the project package importable when the script is run directly.
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fragblast.align import align_pair
from fragblast.segments import find_fragments


def read_fasta(path: Path) -> list[tuple[str, str, str]]:
    """Return ``(identifier, full_header, sequence)`` in file order.

    Tolerates the empty ``>`` line that is present at the end of the supplied
    ``IHSMGC.pep`` file; every real header must have an identifier.
    """
    records: list[tuple[str, str, str]] = []
    current: tuple[str, str] | None = None
    chunks: list[str] = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current is not None:
                    records.append((current[0], current[1], "".join(chunks)))
                fields = line[1:].split()
                if not fields:
                    current = None
                    chunks = []
                    continue
                current = (fields[0], line[1:])
                chunks = []
            elif current is not None:
                chunks.append(line.strip())
            elif line.strip():
                raise ValueError(
                    f"sequence data before first header in {path}"
                )
    if current is not None:
        records.append((current[0], current[1], "".join(chunks)))
    return records


def select_fasta(
    path: Path, wanted: set[str]
) -> list[tuple[str, str, str]]:
    """Single streaming pass: keep only records whose id is in ``wanted``."""
    if not wanted:
        return []
    kept: list[tuple[str, str, str]] = []
    current: tuple[str, str] | None = None
    chunks: list[str] = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current is not None and current[0] in wanted:
                    kept.append((current[0], current[1], "".join(chunks)))
                fields = line[1:].split()
                if not fields:
                    current = None
                    chunks = []
                    continue
                current = (fields[0], line[1:])
                chunks = []
            elif current is not None:
                chunks.append(line.strip())
    if current is not None and current[0] in wanted:
        kept.append((current[0], current[1], "".join(chunks)))
    return kept


def run_diamond(
    *,
    query: Path,
    db_fasta: Path,
    db_dmnd: Path,
    out_tsv: Path,
    sensitivity: str,
    force: bool,
    ident_pct: float,
    cov_pct: float,
    coverage_sequence: str,
    threads: int,
) -> None:
    """Run ``diamond makedb`` + ``blastp`` when outputs do not exist yet."""
    diamond = shutil.which("diamond")
    if diamond is None:
        sys.exit("diamond not found; activate the fragblast conda environment")

    if not db_dmnd.exists() or force:
        print(f"[1/3] diamond makedb -> {db_dmnd}", flush=True)
        subprocess.run(
            [
                diamond,
                "makedb",
                "--in",
                str(db_fasta),
                "--db",
                str(db_dmnd),
                "--threads",
                str(threads),
            ],
            check=True,
        )
    else:
        print(f"[1/3] reuse existing DIAMOND database {db_dmnd}", flush=True)

    if not out_tsv.exists() or force:
        print(f"[2/3] diamond blastp -> {out_tsv}", flush=True)
        outfmt_cols = [
            "qseqid",
            "sseqid",
            "pident",
            "qcovhsp",
            "length",
            "qlen",
            "slen",
            "evalue",
            "bitscore",
        ]
        if coverage_sequence in ("both", "subject"):
            outfmt_cols.append("scovhsp")
        outfmt = "6 " + " ".join(outfmt_cols)
        # Pre-filters are derived from the requested thresholds: a potential
        # fragment with identity > ident and coverage > cov requires at least
        # pident * effective_cov > ident * cov on the DIAMOND HSP.
        pre_ident = max(1, math.floor(ident_pct * cov_pct / 100.0))
        pre_cov = math.floor(cov_pct)
        command = [
            diamond,
            "blastp",
            "--db",
            str(db_dmnd),
            "--query",
            str(query),
            "--out",
            str(out_tsv),
            "--outfmt",
            *outfmt.split(),
            *(
                ["--more-sensitive"]
                if sensitivity == "more-sensitive"
                else []
            ),
            "--id",
            str(pre_ident),
            "--query-cover",
            str(pre_cov),
        ]
        if coverage_sequence in ("both", "subject"):
            command += [
                "--subject-cover",
                str(pre_cov),
            ]
        command += [
            "--threads",
            str(threads),
            "--max-target-seqs",
            "0",  # unlimited targets, so no candidate is hidden by rank
            "--max-hsps",
            "100",
        ]
        subprocess.run(
            command,
            check=True,
        )
    else:
        print(f"[2/3] reuse existing DIAMOND output {out_tsv}", flush=True)


def load_diamond_rows(
    path: Path, coverage_sequence: str
) -> list[dict[str, str | float | int]]:
    rows: list[dict[str, str | float | int]] = []
    need_scov = coverage_sequence in ("both", "subject")
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if need_scov and len(parts) != 10:
                raise ValueError(
                    f"{path} lacks the scovhsp column; re-run DIAMOND with "
                    "coverage_sequence=subject/both (or delete the file and "
                    "use --force-diamond)"
                )
            if len(parts) not in (9, 10):
                raise ValueError(
                    f"unexpected column count in {path}: {len(parts)}"
                )
            q, s, pident, qcov, length, qlen, slen, evalue, bitscore = parts[
                :9
            ]
            row: dict[str, str | float | int] = {
                "query": q,
                "subject": s,
                "pident": float(pident),
                "qcov": float(qcov),
                "length": int(length),
                "qlen": int(qlen),
                "slen": int(slen),
                "evalue": float(evalue),
                "bitscore": float(bitscore),
            }
            if len(parts) == 10:
                row["scov"] = float(parts[9])
            rows.append(row)
    return rows


def coverage_value(
    row: dict[str, str | float | int], coverage_sequence: str
) -> float:
    """Coverage used for filtering, following --coverage_sequence."""
    if coverage_sequence == "query":
        return float(row["qcov"])
    if coverage_sequence == "subject":
        return float(row.get("scov", 0.0))
    return min(float(row["qcov"]), float(row.get("scov", 0.0)))


def coverage_ok(
    row: dict[str, str | float | int],
    coverage_sequence: str,
    cov_pct: float,
) -> bool:
    """True when every sequence required by coverage_sequence passes ``cov``."""
    if coverage_sequence in ("query", "both"):
        if float(row["qcov"]) <= cov_pct:
            return False
    if coverage_sequence in ("subject", "both"):
        if float(row.get("scov", 0.0)) <= cov_pct:
            return False
    return True


def write_fasta(path: Path, records: list[tuple[str, str, str]], ids) -> None:
    by_id = {rec[0]: rec for rec in records}
    with open(path, "w", encoding="utf-8") as handle:
        for seq_id in ids:
            _, header, seq = by_id[seq_id]
            handle.write(f">{header}\n")
            for pos in range(0, len(seq), 80):
                handle.write(seq[pos : pos + 80] + "\n")


def verify_with_fragblast(
    *,
    query_records: list[tuple[str, str, str]],
    db_records: list[tuple[str, str, str]],
    candidate_pairs: list[dict[str, str | float | int]],
    ident_pct: float,
    qcov_pct: float,
    scov_pct: float,
    out_tsv: Path,
    workers: int = 1,
) -> set[str]:
    """Re-align candidate pairs with parasail SW and enumerate fragments."""
    query_by_id = {rec[0]: rec for rec in query_records}
    subject_by_id = {rec[0]: rec for rec in db_records}
    confirmed: set[str] = set()
    tasks = [
        (
            row["query"],
            row["subject"],
            query_by_id[row["query"]][2],
            subject_by_id[row["subject"]][2],
            ident_pct,
            qcov_pct,
            scov_pct,
        )
        for row in candidate_pairs
    ]
    header = (
        "query_id\tsubject_id\tq_start\tq_end\ts_start\ts_end"
        "\taligned_len\tidentity\tqcov\tscov\n"
    )
    with open(out_tsv, "w", encoding="utf-8") as handle:
        handle.write(header)
        if workers > 1 and len(tasks) > 1:
            ctx = mp.get_context("spawn")
            with ctx.Pool(workers) as pool:
                for query_id, fragments in pool.imap_unordered(
                    _verify_one_pair, tasks, chunksize=16
                ):
                    for frag in fragments:
                        handle.write(frag)
                    if fragments:
                        confirmed.add(query_id)
        else:
            for task in tasks:
                query_id, fragments = _verify_one_pair(task)
                for frag in fragments:
                    handle.write(frag)
                if fragments:
                    confirmed.add(query_id)
    return confirmed


def _verify_one_pair(task: tuple) -> tuple[str, list[str]]:
    """Worker for verify_with_fragblast (returns one TSV row per fragment)."""
    (
        query_id,
        subject_id,
        query_seq,
        subject_seq,
        ident_pct,
        qcov_pct,
        scov_pct,
    ) = task
    alignment = align_pair(query_seq, subject_seq)
    if alignment is None:
        return query_id, []
    fragments = find_fragments(
        alignment.columns,
        len(query_seq),
        len(subject_seq),
        identity_pct=ident_pct,
        qcov_pct=qcov_pct,
        scov_pct=scov_pct,
        q_offset=alignment.q_offset,
        s_offset=alignment.s_offset,
    )
    rows: list[str] = []
    for frag in fragments:
        rows.append(
            "\t".join(
                str(v)
                for v in (
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
                )
            )
            + "\n"
        )
    return query_id, rows


def write_summary(
    path: Path,
    *,
    n_diamond_rows: int,
    solid_rows: list,
    solid_ids: list[str],
    potential_ids: list[str],
    candidate_rows: list,
    confirmed_ids: list[str],
) -> None:
    solid_subjects = sorted({r["subject"] for r in solid_rows})
    potential_subjects = sorted({r["subject"] for r in candidate_rows})
    lines = [
        "# DIAMOND + FragBLAST: solid/potential/new ARG summary",
        f"diamond_hsp_rows={n_diamond_rows}",
        f"solid_rows={len(solid_rows)}",
        f"solid_query_sequences={len(solid_ids)}",
        f"solid_subject_sequences={len(solid_subjects)}",
        f"potential_query_sequences={len(potential_ids)}",
        f"candidate_pairs_aligned={len(candidate_rows)}",
        f"candidate_subject_sequences={len(potential_subjects)}",
        f"fragblast_confirmed_queries={len(confirmed_ids)}",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines), flush=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--query", type=Path, default=ROOT / "data" / "IHSMGC.pep")
    parser.add_argument("--db", type=Path, default=ROOT / "data" / "SARG.fasta")
    parser.add_argument("--outdir", type=Path, default=ROOT / "results")
    parser.add_argument("--ident", type=float, default=80.0)
    parser.add_argument("--cov", type=float, default=80.0)
    parser.add_argument(
        "--coverage_sequence",
        choices=("query", "subject", "both"),
        default="query",
        help=(
            "which sequence(s) the --cov threshold applies to: query, "
            "subject or both (default: query)"
        ),
    )
    parser.add_argument(
        "--sensitivity",
        choices=("default", "more-sensitive"),
        default="more-sensitive",
    )
    parser.add_argument(
        "--force-diamond",
        action="store_true",
        help="re-run diamond even if its outputs already exist",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="CPU threads for DIAMOND (default: 4)",
    )
    parser.add_argument(
        "--max-subjects-per-query",
        type=int,
        default=3,
        help=(
            "keep at most this many best DIAMOND subjects per potential "
            "query for FragBLAST (default: 3)"
        ),
    )
    parser.add_argument(
        "--fragblast-threads",
        type=int,
        default=4,
        help="parallel workers for FragBLAST verification (default: 4)",
    )
    args = parser.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)
    query_name = args.query.name
    for suffix in (
        ".pep.gz",
        ".fa.gz",
        ".fasta.gz",
        ".faa.gz",
        ".pep",
        ".fa",
        ".fasta",
        ".faa",
    ):
        if query_name.endswith(suffix):
            query_name = query_name[: -len(suffix)]
            break
    query_base = query_name
    db_dmnd = args.outdir / f"{args.db.stem}.dmnd"
    diamond_tsv = args.outdir / f"{query_base}_vs_{args.db.stem}.tsv"
    solid_fasta = args.outdir / "solid_ARGs.fasta"
    potential_fasta = args.outdir / "potential_ARGs.fasta"
    potential_pairs_tsv = args.outdir / "potential_pairs.tsv"
    fragblast_tsv = args.outdir / "fragblast_potential.tsv"
    new_fasta = args.outdir / "new_ARGs.fasta"
    summary = args.outdir / "arg_pipeline_summary.txt"

    run_diamond(
        query=args.query,
        db_fasta=args.db,
        db_dmnd=db_dmnd,
        out_tsv=diamond_tsv,
        sensitivity=args.sensitivity,
        force=args.force_diamond,
        ident_pct=args.ident,
        cov_pct=args.cov,
        coverage_sequence=args.coverage_sequence,
        threads=args.threads,
    )

    print("[3/3] filtering and FragBLAST verification", flush=True)
    db_records = read_fasta(args.db)
    if not db_records:
        parser.error(f"no sequences found in db file {args.db}")

    rows = load_diamond_rows(diamond_tsv, args.coverage_sequence)
    solid_rows = [
        r
        for r in rows
        if r["pident"] > args.ident
        and coverage_ok(r, args.coverage_sequence, args.cov)
    ]

    def product(r) -> float:
        return (r["pident"] / 100.0) * (
            coverage_value(r, args.coverage_sequence) / 100.0
        )

    solid_query_set = {r["query"] for r in solid_rows}
    # A candidate must already cover >cov of the relevant sequence(s): no
    # embedded fragment can exceed the coverage of the whole HSP. It must
    # also satisfy pident * coverage > ident * cov, a necessary condition for
    # a qualifying fragment to exist inside the HSP.
    product_threshold = (args.ident / 100.0) * (args.cov / 100.0)
    product_rows_all = [
        r
        for r in rows
        if r["query"] not in solid_query_set
        and coverage_ok(r, args.coverage_sequence, args.cov)
        and product(r) > product_threshold
    ]
    product_query_set = {r["query"] for r in product_rows_all}
    potential_query_set = product_query_set - solid_query_set

    # One streaming pass through the (possibly huge) query FASTA keeps only
    # the records that will be written or re-aligned.
    selected_records = select_fasta(
        args.query, solid_query_set | potential_query_set
    )
    solid_records = [r for r in selected_records if r[0] in solid_query_set]
    potential_records = [
        r for r in selected_records if r[0] in potential_query_set
    ]
    solid_query_ids = [r[0] for r in solid_records]
    potential_query_ids = [r[0] for r in potential_records]
    # One local alignment per query-subject pair is enough: if any DIAMOND
    # HSP met the product cutoff we keep only the best HSP for that pair.
    pair_best: dict[tuple[str, str], dict] = {}
    found_query_set = {r[0] for r in potential_records}
    for r in product_rows_all:
        if r["query"] not in found_query_set:
            continue
        key = (r["query"], r["subject"])
        if key not in pair_best or product(r) > product(pair_best[key]):
            pair_best[key] = r
    by_query: dict[str, list[dict]] = {}
    for key, best_row in pair_best.items():
        by_query.setdefault(key[0], []).append(best_row)
    candidate_rows: list[dict] = []
    for rows_for_query in by_query.values():
        rows_for_query.sort(key=product, reverse=True)
        candidate_rows.extend(
            rows_for_query[: args.max_subjects_per_query]
        )

    write_fasta(solid_fasta, solid_records, solid_query_ids)
    write_fasta(potential_fasta, potential_records, potential_query_ids)
    with open(potential_pairs_tsv, "w", encoding="utf-8") as handle:
        cols = [
            "query_id",
            "subject_id",
            "pident",
            "qcovhsp",
            "effective_cov",
            "product",
            "length",
            "qlen",
            "slen",
            "evalue",
            "bitscore",
        ]
        if args.coverage_sequence in ("both", "subject"):
            cols.insert(4, "scovhsp")
        handle.write("\t".join(cols) + "\n")
        for r in candidate_rows:
            values = [
                r["query"],
                r["subject"],
                r["pident"],
                r["qcov"],
            ]
            if args.coverage_sequence in ("both", "subject"):
                values.append(r.get("scov", ""))
            values.extend(
                [
                    coverage_value(r, args.coverage_sequence),
                    product(r),
                    r["length"],
                    r["qlen"],
                    r["slen"],
                    r["evalue"],
                    r["bitscore"],
                ]
            )
            handle.write(
                "\t".join(
                    str(v)
                    for v in values
                )
                + "\n"
            )

    qcov_pct = args.cov if args.coverage_sequence in ("query", "both") else 0.0
    scov_pct = args.cov if args.coverage_sequence in ("subject", "both") else 0.0
    confirmed = verify_with_fragblast(
        query_records=potential_records,
        db_records=db_records,
        candidate_pairs=candidate_rows,
        ident_pct=args.ident,
        qcov_pct=qcov_pct,
        scov_pct=scov_pct,
        out_tsv=fragblast_tsv,
        workers=args.fragblast_threads,
    )
    confirmed_records = [
        r for r in potential_records if r[0] in confirmed
    ]
    confirmed_ids = [r[0] for r in confirmed_records]
    write_fasta(new_fasta, confirmed_records, confirmed_ids)

    write_summary(
        summary,
        n_diamond_rows=len(rows),
        solid_rows=solid_rows,
        solid_ids=solid_query_ids,
        potential_ids=potential_query_ids,
        candidate_rows=candidate_rows,
        confirmed_ids=confirmed_ids,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
