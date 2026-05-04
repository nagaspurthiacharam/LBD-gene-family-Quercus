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
│
├── README.md
├── LICENSE
├── 01_hmmer_identification.sh
│
├── 02_gene_expression/
│   ├── README.md
│   ├── DEG_analysis.R
│   ├── extract_QrLBD_counts.R
│   ├── Heatmap_zscore_QrLBD.R
│   ├── Heatmap_zscore_QsLBD.R
│   └── PE_paired_end/
│       ├── 0_downloadSRA_PE.sh
│       ├── 1_QualityCheck_PE.sh
│       ├── 2_Trimmomatic_PE.sh
│       └── 3_hisat2_mapping_PE.sh
│
├── 03_collinearity/
│   ├── collinearity_workflow.sh
│   └── plot_synteny_all.py
│
├── 04_cis_promoter/
│   ├── Cis_tile_heatmap.R
│   ├── Cis_curated_QrLBD42.xlsx
│   └── Cis_curated_QsLBD37.xlsx
│
└── data/
    ├── LBD_gene_IDs.txt
    ├── QrLBD_gene_catalogue.csv          
    ├── QrLBD_phylogeny_order.txt
    ├── QsLBD_phylogeny_order.txt
    ├── QrLBD_proteins.fasta              (42 sequences)
    ├── QsLBD_proteins.fasta              (37 sequences)
    ├── QrLBD_promoters_2kb.fasta
    └── QsLBD_promoters_2kb.fasta
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

PlantCARE prediction across the 2-kb upstream region of every LBD locus, curated to 22 informative motifs grouped under five functional categories (hormone-, stress-, light/circadian-, promoter-core-, transcription-factor-binding-related). The R script `Cis_tile_heatmap.R` produces:

- **Main figure** (`Cis_main_summary.png`): five-category totals for both species side-by-side with a shared colour scale for direct interspecific comparison.
- **Supplementary figures** (`Cis_QrLBD_tile.png`, `Cis_QsLBD_tile.png`): per-element detailed matrices, one per species.

Tile-heatmap convention (raw counts + heat-colour fill) follows the *Brassica napus* PSK family analysis (Zhang et al., 2026) and the LBD genome-wide surveys in poplar, bamboo, and grapevine.

## Genome Resources

| Species | Accession | Source |
|---|---|---|
| *Quercus rubra* v2.1 | — | Phytozome v13 |
| *Quercus suber* | GCF_002906115.3 | NCBI RefSeq |
| *Populus trichocarpa* v4.1 | GCF_000002775.5 | NCBI RefSeq |
| *Arabidopsis thaliana* TAIR10.1 | GCA_000001735.2 | NCBI GenBank |
| *Oryza sativa* indica 93-11 | GCA_000004655.2 | NCBI GenBank |
| *Populus alba* | — | Phytozome v13 |

## Software & Citations

| Tool | Version | Citation |
|---|---|---|
| HMMER | 3.3.2 | Eddy (2011) PLoS Comput Biol 7:e1002195 |
| NCBI BLAST+ | 2.12 | Camacho et al. (2009) BMC Bioinformatics 10:421 |
| NCBI CDD / CD-Search | — | Marchler-Bauer et al. (2017) Nucleic Acids Res 45:D200 |
| HISAT2 | 2.2.1 | Kim et al. (2019) Nat Biotechnol 37:907 |
| StringTie | 2.2.1 | Pertea et al. (2015) Nat Biotechnol 33:290 |
| Trimmomatic | 0.39 | Bolger et al. (2014) Bioinformatics 30:2114 |
| FastQC | 0.12 | Andrews (2010) Babraham Bioinformatics |
| pheatmap | 1.0.12 | Kolde (2019) CRAN |
| ggplot2 | 3.4 | Wickham (2016) Springer |
| patchwork | 1.1 | Pedersen (2024) CRAN |
| MCScanX | — | Wang et al. (2012) Nucleic Acids Res 40:e49 |
| matplotlib | 3.7 | Hunter (2007) Comput Sci Eng 9:90 |
| MEGA 7 | 7.0 | Kumar et al. (2016) Mol Biol Evol 33:1870 |
| MEME Suite | 5.5 | Bailey et al. (2015) Nucleic Acids Res 43:W39 |
| TBtools-II | 2.0 | Chen et al. (2023) Mol Plant 16:1733 |
| iTOL | 7 | Letunic & Bork (2024) Nucleic Acids Res 52:W78 |
| Blast2GO | 6.0 | Conesa et al. (2005) Bioinformatics 21:3674 |
| STRING | 11.5 | Szklarczyk et al. (2021) Nucleic Acids Res 49:D605 |
| PlantCARE | — | Lescot et al. (2002) Nucleic Acids Res 30:325 |
| MG2C | v2.1 | Chao et al. (2015, 2021) |

## Requirements

```bash
sudo apt-get install ncbi-blast+ hisat2 stringtie trimmomatic fastqc samtools
conda install -c bioconda hmmer
git clone https://github.com/wyp1125/MCScanX.git && cd MCScanX && make
Rscript -e 'install.packages(c("pheatmap","RColorBrewer","ggplot2","readxl","dplyr","tidyr","scales","patchwork"))'
pip install matplotlib pandas openpyxl
```

## Setup

Replace `<YOUR_HPC_ID>` in shell scripts with your HPC username before running.

## Acknowledgments

The RNA-seq analysis scripts in `02_gene_expression/` were adapted from the Functional Genomics course (BIOL7180) pipeline originally developed by **Dr. Rita Graze**, Auburn University.

## License

MIT License. Please cite the tools listed above when using this pipeline in publications.
