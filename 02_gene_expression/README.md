# 02_gene_expression — RNA-seq pipeline & LBD expression heatmaps

This folder runs the complete RNA-seq analysis used in the manuscript: from raw SRA reads, through QC and mapping, all the way to the LBD-family z-score expression heatmaps.

## Datasets

| Dataset | Species | BioProject | Library | Tissues | Samples |
|---|---|---|---|---|---|
| *Q. rubra* multitissue | *Q. rubra* | PRJNA273270 | Paired-end | emerging leaf/bud, immature twig, twig late-growth, mature leaf, dormant twig, developing acorn (two parental genotypes SM1/SM2) | 11 |
| *Q. suber* multitissue | *Q. suber* | PRJNA392919 | Paired-end | phellem, xylem, leaf, inner bark, pollen | 5 |

## Folder layout

```
02_gene_expression/
├── README.md                         ← this file
├── PE_paired_end/                    ← upstream HPC pipeline (raw reads → counts)
│   ├── 0_downloadSRA_PE.sh
│   ├── 1_QualityCheck_PE.sh
│   ├── 2_Trimmomatic_PE.sh
│   └── 3_hisat2_mapping_PE.sh
├── DEG_analysis.R                    ← generic DESeq2 reference workflow
├── extract_QrLBD_counts.R            ← regenerate ../data/QrLBD_counts_42_primary.csv
├── extract_QsLBD_counts.R            ← regenerate ../data/QsLBD_counts_37_primary.csv
├── Heatmap_zscore_QrLBD.R            ← Figure 4 (Q. rubra)
└── Heatmap_zscore_QsLBD.R            ← Figure 5 (Q. suber)
```

## Pipeline overview

```
SRA reads (PRJNA273270 / PRJNA392919)
        │
        ▼   PE_paired_end/ shell scripts on the HPC
   FastQC ──► Trimmomatic ──► HISAT2 ──► StringTie ──► prepDE.py3
        │
        ▼
   gene_count_matrix.csv         (whole-genome, NOT shipped)
        │
        ▼   extract_Qr/QsLBD_counts.R
   ../data/QrLBD_counts_42_primary.csv       ← shipped
   ../data/QsLBD_counts_37_primary.csv       ← shipped
        │
        ▼   Heatmap_zscore_Qr/QsLBD.R
   Figures 4 & 5 (z-score expression heatmaps)
```

## 1. Upstream HPC pipeline — `PE_paired_end/`

Four numbered shell scripts that take raw SRA accessions and produce a per-sample BAM ready for StringTie. Replace `<YOUR_HPC_ID>` with your HPC username before running.

| # | Script | Purpose |
|---|---|---|
| 0 | `0_downloadSRA_PE.sh` | `prefetch` + `fastq-dump --split-files` for every SRR accession |
| 1 | `1_QualityCheck_PE.sh` | FastQC pre- and post-trim |
| 2 | `2_Trimmomatic_PE.sh` | adapter / quality trimming, PE mode (4 output files per sample) |
| 3 | `3_hisat2_mapping_PE.sh` | HISAT2 splice-aware alignment, sort to coord-sorted BAM |

After mapping, StringTie is run per-sample (`-eB` mode against the reference annotation), and `prepDE.py3` collapses the per-sample `.ctab` outputs into a single whole-genome `gene_count_matrix.csv` (33,247 genes for *Q. rubra* / 43,163 for *Q. suber*). This whole-genome matrix is the input to the LBD subset extractors below; it is **not** shipped with the repository (it is regenerable from public SRA reads).

## 2. LBD subset extractors

Both extractors read a local whole-genome `gene_count_matrix.csv` (placed next to the script) and write the LBD-family-only subset into `../data/`. On a fresh checkout — i.e. if no whole-genome matrix is present — the scripts print an informative message and exit cleanly: the LBD-only subsets are already shipped in `../data/`, so the heatmap step works straight out of the box.

| Script | Reads | Writes | Rows |
|---|---|---|---|
| `extract_QrLBD_counts.R` | `../data/QrLBD_gene_catalogue.csv` (catalogue) + local `gene_count_matrix.csv` | `../data/QrLBD_counts_42_primary.csv` | 42 |
| `extract_QsLBD_counts.R` | embedded 37-row QsLBD ↔ NCBI LOC mapping + local `gene_count_matrix.csv` | `../data/QsLBD_counts_37_primary.csv` | 37 |

The Qr extractor uses the catalogue's `Sub-isoform of` column to keep primaries only (42 of 47 rows). The Qs extractor uses an embedded QsLBD1–QsLBD37 ↔ `gene-LOCxxxxxxxxx|LOCxxxxxxxxx` mapping table (the analogue of the Qr catalogue).

Run from this folder:

```bash
Rscript extract_QrLBD_counts.R
Rscript extract_QsLBD_counts.R
```

## 3. Z-score heatmaps — `Heatmap_zscore_*.R`

Both scripts read the shipped LBD-only subset in `../data/` (so they run end-to-end with no extra downloads), apply identical statistics, and render the figure via `pheatmap`.

```r
# Q. rubra
COUNTS <- "../data/QrLBD_counts_42_primary.csv"
ORDER  <- "../data/QrLBD_phylogeny_order.txt"

# Q. suber
COUNTS <- "../data/QsLBD_counts_37_primary.csv"
ORDER  <- "../data/QsLBD_phylogeny_order.txt"
```

Common statistics in both heatmap scripts:

1. DESeq2-style prefilter — drop genes with `rowSums(counts) <= 20` across the full library set (matches the manuscript Methods).
2. Median-of-ratios library-size normalization (DESeq2 `estimateSizeFactors` semantics).
3. Average across replicates per tissue — *Q. rubra* only; *Q. suber* has one library per tissue.
4. Per-gene z-score, capped at ±2 to compress outliers without truncating the colour scale.
5. Row order driven by the mid-point-rooted phylogeny leaf order (`Qr/QsLBD_phylogeny_order.txt`), so the heatmap, the gene-structure / motif diagram, and the cis-element heatmap all share the same biological row order.

Outputs: `Heatmap_zscore_QrLBD.png` (Figure 4) and `Heatmap_zscore_QsLBD.png` (Figure 5).

## 4. `DEG_analysis.R` — generic DESeq2 reference

A standalone reference workflow showing the canonical DESeq2 differential-expression approach (design formula, `DESeq()`, `results()`, MA-plot, volcano). It is **not** invoked by the heatmap pipeline; it is included so reviewers can see the DESeq2 conventions that the heatmap scripts borrow (prefilter rule, median-of-ratios sizing).

## Sample columns in the shipped subsets

`QrLBD_counts_42_primary.csv` (42 rows × 14 cols):

```
QrLBD, rna_id, gene_id,
SRR1778917, SRR1778919, SRR1778921, SRR1778922, SRR1778923,
SRR1778925, SRR1778927, SRR1778929, SRR1778930, SRR1778931, SRR1778932
```

`QsLBD_counts_37_primary.csv` (37 rows × 7 cols):

```
QsLBD, gene_id,
SRR5986737, SRR5986738, SRR5986739, SRR5986740, SRR5986741
```

## Software

| Tool | Version | Citation |
|---|---|---|
| FastQC | 0.12 | Andrews (2010) Babraham Bioinformatics |
| Trimmomatic | 0.39 | Bolger et al. (2014) Bioinformatics 30:2114 |
| HISAT2 | 2.2.1 | Kim et al. (2019) Nat Biotechnol 37:907 |
| StringTie | 2.2.1 | Pertea et al. (2015) Nat Biotechnol 33:290 |
| prepDE.py3 | (StringTie auxiliary) | shipped with StringTie |
| R / DESeq2 | 4.x / 1.40 | Love et al. (2014) Genome Biol 15:550 |
| pheatmap | 1.0.12 | Kolde (2019) CRAN |

## Setup

Replace `<YOUR_HPC_ID>` in `PE_paired_end/*.sh` with your HPC username before running. The R scripts have no HPC dependency and run on any machine with R + DESeq2 + pheatmap installed.

## Acknowledgment

The shell scripts in `PE_paired_end/` were adapted from the Functional Genomics course (BIOL7180) pipeline originally developed by **Dr. Rita Graze**, Auburn University.
