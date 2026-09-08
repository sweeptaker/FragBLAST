# FragBLAST

FragBLAST is a fragment-level protein alignment workflow for detecting
antibiotic resistance genes (ARGs). It combines DIAMOND for fast candidate
retrieval with a local re-alignment step (parasail / Smith–Waterman) that
enumerates maximal high-identity fragments along a single alignment path.

## Why FragBLAST?

DIAMOND and BLAST report whole HSPs, and hard thresholds such as
"identity > 80% and query coverage > 80%" are applied to averages over each
entire HSP. A long HSP can therefore pass or fail a threshold even though it
contains an internal segment with much higher identity.

Example: a query-subject pair may have an overall HSP of 60% identity and 60%
query coverage, while a fragment inside that HSP has > 80% identity and > 80%
query coverage. Standard hard-threshold filtering reports only the average
properties of the whole HSP and never exposes the embedded high-identity
fragment.

FragBLAST addresses this by:

1. using DIAMOND to find candidate query-subject pairs;
2. re-aligning each candidate pair with full local alignment (BLOSUM62,
   affine gaps);
3. enumerating maximal contiguous fragments of the alignment path whose
   identity, query coverage and subject coverage exceed user-defined
   thresholds;
4. reporting those fragments as additional candidate ARG regions.

The risk of underestimating embedded high-identity fragments is higher when a
lenient threshold is used, because long HSPs with mixed identity are more
likely to be accepted as whole alignments, hiding the high-identity segment
inside them.

## Datasets and ARG results

The workflow was applied to three protein/gene catalogs:

| Dataset | Source | Sample type | DIAMOND threshold | FragBLAST threshold |
|---|---|---|---|---|
| IHSMGC | Human skin microbiome gene catalog | Human | identity > 80%, qcov > 80% | identity > 80%, qcov > 80% |
| IGC | Integrated human gut microbiome gene catalog | Human | identity > 80%, qcov > 80% | identity > 80%, qcov > 80% |
| Tara Ocean | Tara Oceans expedition proteins | Environment | identity > 70%, qcov > 40% | identity > 70%, qcov > 40% |

The two human catalogs were screened with the strict 80%/80% rule. The
environmental Tara Ocean catalog was screened with a more lenient 70%/40% rule
because of the larger diversity and lower completeness of environmental
sequence data. **Results obtained with different thresholds should not be
compared directly.**

### ARG counts

For each dataset, "Existing ARGs" are the DIAMOND hits that already pass the
hard threshold, and "New ARGs" are the additional sequences confirmed by
FragBLAST from the remaining candidate pool.

| Dataset | Threshold | Existing ARGs | New ARGs |
|---|---|---:|---:|
| IHSMGC (human skin) | 80% / 80% | 5,534 | 857 |
| IGC (human gut) | 80% / 80% | 1,176 | 138 |
| Tara Ocean (environment) | 70% / 40% | 1,157 | 3,375 |
| Tara Ocean (environment, re-analyzed) | 80% / 80% | 118 | 96 |

When Tara Ocean is re-analyzed with the same 80%/80% rule as the human
catalogs, the number of detected ARGs drops sharply, illustrating how strongly
the threshold choice affects both detection and interpretation.

## Repository layout

```text
.
├── fragblast/                  # core fragment enumeration library
│   ├── align.py                # parasail local alignment wrapper
│   ├── segments.py             # maximal fragment enumeration
│   └── cli.py                  # command-line interface
├── scripts/
│   ├── arg_pipeline.py         # DIAMOND + FragBLAST ARG workflow
│   └── plot_arg_comparison_v2.R  # ggplot2 summary figures
├── tests/                      # pytest suite
├── figures/                    # English summary figures
└── environment.yml             # conda environment
```

## Installation

The recommended setup uses conda:

```bash
conda env create -f environment.yml
conda activate fragblast
```

Install the parasail Python bindings. The 1.3.4 source distribution requires a
`glibtoolize` alias on macOS:

```bash
mkdir -p .build-bin
ln -sf "$(command -v libtoolize)" .build-bin/glibtoolize
ln -sf "$(command -v libtool)" .build-bin/glibtool
PATH="$PWD/.build-bin:$PATH" python -m pip install parasail==1.3.4
```

Install DIAMOND if it is not already available:

```bash
conda install -c bioconda diamond
```

Run the test suite:

```bash
python -m pytest
```

## Usage

### FragBLAST command line

```bash
python -m fragblast \
  --query query.fa \
  --db database.fa \
  --ident 80 \
  --cov 80 \
  --coverage_sequence both \
  --out hits.tsv
```

Output columns:

`query_id  subject_id  q_start  q_end  s_start  s_end  aligned_len  identity  qcov  scov`

### End-to-end ARG pipeline

```bash
python scripts/arg_pipeline.py \
  --query data/query.fa \
  --db data/SARG.fasta \
  --outdir results \
  --ident 80 \
  --cov 80 \
  --product 0.64 \
  --diamond-id 64 \
  --diamond-qcov 64 \
  --threads 4 \
  --fragblast-threads 4 \
  --max-subjects-per-query 3
```

The pipeline writes:

- `solid_ARGs.fasta`: DIAMOND hits that pass the hard threshold;
- `potential_ARGs.fasta`: remaining candidates that pass the relaxed product
  threshold;
- `new_ARGs.fasta`: candidates additionally confirmed by FragBLAST;
- `arg_pipeline_summary.txt`: sequence counts and summary statistics.

## Figures

![ARG counts under the 80/80 rule](figures/arg_counts_80_80.png)

![Effect of the detection criteria on Tara Ocean ARGs](figures/tara_standard_comparison.png)

## References

- Li, Z. et al. Characterization of the human skin resistome and
  identification of two microbiota cutotypes. Microbiome 9, 47 (2021).
- Li, J. et al. An integrated catalog of reference genes in the human gut
  microbiome. Nature Biotechnology 32, 834–841 (2014).
- Delmont, T. O. et al. Nitrogen-fixing populations of Planctomycetes and
  Proteobacteria are abundant in surface ocean metagenomes. Nature
  Microbiology 3, 804–813 (2018).
