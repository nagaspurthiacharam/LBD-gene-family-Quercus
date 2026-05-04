# Genome-Wide Identification and Characterization of the LBD Gene Family in *Quercus rubra* and *Quercus suber*

**Author:** Naga Spurthi Acharam
**Affiliation:** College of Forestry, Wildlife and Environment, Auburn University
**Date:** 2026

---

## Overview

This repository contains the bioinformatics pipeline for the genome-wide identification, phylogenetic classification, gene-expression profiling, synteny analysis, and cis-regulatory element characterization of the Lateral Organ Boundary Domain (LBD) gene family in two oak species with contrasting crown architectures:

- ***Quercus rubra*** (Northern Red Oak) - excurrent growth, strong apical dominance
- ***Quercus suber*** (Cork Oak) - decurrent growth, weak apical dominance, thick cork-producing phellem

Final analytical catalogue: **42 QrLBD primary loci + 37 QsLBD primary loci = 79 genes total** (five sub-isoforms tabulated separately for transparency but excluded from downstream analyses).

## Repository Structure

```
LBD-gene-family-Quercus/
|
|-- README.md
|-- 01_hmmer_identification.sh         HMMER + BLASTP + NCBI CDD gene identification
|-- 02_gene_expression/
|   |-- README.md                      Pipeline overview (paired-end)
|   |-- PE_paired_end/                 SRA download, QC, mapping pipeline
|   |   |-- 0_downloadSRA_PE.sh
|   |   |-- 1_QualityCheck_PE.sh
|   |   |-- 2_Trimmomatic_PE.sh
|   |   |-- 3_hisat2_mapping_PE.sh
|   |-- extract_QrLBD_counts.R         Pull 42-gene counts from full matrix
|   |-- Heatmap_zscore_QrLBD.R         z-score heatmap, Q. rubra (PRJNA273270)
|   |-- Heatmap_zscore_QsLBD.R         z-score heatmap, Q. suber (PRJNA392919)
|-- 03_collinearity/
|   |-- collinearity_workflow.sh       BLASTP + MCScanX pipeline
|   |-- plot_synteny_all.py            Unified four-comparison synteny plotter
|-- 04_cis_promoter/
|   |-- Cis_tile_heatmap.R             Tile heatmap of cis-elements (main + supp)
|   |-- Cis_curated_QrLBD42.xlsx       Curated per-element counts, Q. rubra
|   |-- Cis_curated_QsLBD37.xlsx       Curated per-element counts, Q. suber
|-- data/
    |-- LBD_highlight_QrQs.txt         79 LBD gene IDs (42 Qr + 37 Qs primary)
    |-- QrLBD_master.csv               Canonical 42-gene catalogue + 5 sub-isoforms
    |-- QrLBD_phylo_order.txt          Mid-point-rooted Qr leaf order, 42 entries
    |-- QsLBD_phylo_order.txt          Mid-point-rooted Qs leaf order, 37 entries
    |-- QrLBD_42_proteins.fasta        Q. rubra LBD proteome
    |-- QsLBD_37_proteins.fasta        Q. suber LBD proteome
```

## Analysis Pipeline

| Step | Analysis | Tools | Scripts |
|------|----------|-------|---------|
| 1 | LBD gene identification | HMMER 3.3 (PF03195), BLASTP, NCBI CDD, Pfam web, SMART | `01_hmmer_identification.sh` |
| 2 | Gene expression profiling | HISAT2, StringTie, prepDE.py3, R | `02_gene_expression/` |
| 3 | Genome-wide collinearity | BLASTP, MCScanX, matplotlib | `03_collinearity/` |
| 4 | Cis-regulatory elements | PlantCARE, ggplot2 (tile heatmap) | `04_cis_promoter/` |

## 1. Gene Identification (`01_hmmer_identification.sh`)

Three-step identification pipeline: (i) HMMER v3.3.2 `hmmsearch` of the LOB-domain HMM profile (Pfam accession PF03195.20) against each oak proteome at an E-value threshold of 1 x 10^-5; (ii) reciprocal BLASTP v2.12 search using all 43 *A. thaliana* LBD proteins as queries against each oak proteome (E-value 1 x 10^-5; soft-masking; 250 max hits per query) to catch any divergent Class II members missed by the HMM; (iii) verification of every merged candidate through the NCBI Conserved Domain Database via the Batch Web CD-Search Tool with default parameters - sequences not returning a positive PF03195 (LOB) hit, or retaining only a partial C-block without GAS or coiled-coil signatures, were discarded. All retained sequences cross-validated through the Pfam web server and SMART.

**Results:** 42 *Q. rubra* primaries (QrLBD1-QrLBD42), 37 *Q. suber* primaries (QsLBD1-QsLBD37); five sub-isoforms (QrLBD1.1, QrLBD19.1, QrLBD19.2, QrLBD22.1, QrLBD40.1) flagged at the ProtParam quality-control step and excluded from downstream analyses.

## 2. Gene Expression (`02_gene_expression/`)

RNA-seq analysis of two publicly available paired-end datasets:

- ***Q. rubra* multitissue** (PRJNA273270, six tissues: emerging leaf/bud, immature twig, twig late-growth, mature leaf, dormant twig, developing acorn; two parental genotypes Mother SM1 / Father SM2)
- ***Q. suber* multitissue** (PRJNA392919, five tissues: phellem, xylem, leaf, inner bark, pollen)

**Pipeline:** FastQC -> Trimmomatic -> HISAT2 -> StringTie -> prepDE.py3 -> R

**R analysis (per species):**
1. `extract_QrLBD_counts.R` joins `gene_count_matrix.csv` with `qrlbd_to_rgq29_mapping_NEW.csv` (rna_id <-> gene_id <-> QrLBDs_new) on gene_id, drops sub-isoform rows, writes a 42-gene primary CSV.
2. `Heatmap_zscore_QrLBD.R` and `Heatmap_zscore_QsLBD.R` apply DESeq2-style prefiltering (`rowSums(counts) > 20` on the full matrix), median-of-ratios library-size normalization, average across replicates per tissue (Q. rubra only; Q. suber has one library per tissue), per-gene z-score capped at +/- 2, render `pheatmap` with rows ordered by mid-point-rooted phylogeny.

## 3. Collinearity Analysis (`03_collinearity/`)

Genome-wide synteny detection across four species comparisons via MCScanX. The unified `plot_synteny_all.py` reads the four `.gff` + `.collinearity` outputs and produces dual-track ribbon plots with grey ribbons for all syntenic pairs and red ribbons for LBD-anchored pairs (using `data/LBD_highlight_QrQs.txt` for highlight membership).

| Comparison | Collinear blocks | Gene pairs | LBD pairs |
|---|---|---|---|
| *Q. rubra* vs *Q. suber* | 1,057 | 11,012 | 14 |
| *Q. rubra* vs *P. trichocarpa* | 850 | 19,413 | 40 |
| *Q. rubra* vs *A. thaliana* | 711 | 8,916 | 14 |
| *Q. rubra* vs *O. sativa* | 319 | 2,833 | 9 |

## 4. Cis-Regulatory Elements (`04_cis_promoter/`)

PlantCARE prediction across the 2-kb upstream region of every LBD locus, curated to 22 informative motifs grouped under five functional categories (hormone-, stress-, light/circadian-, promoter-core-, transcription-factor-binding-related). The R scr