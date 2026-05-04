# Gene Expression Analysis

RNA-seq analysis pipeline for profiling LBD gene expression across tissues.

## Datasets

| Dataset | Species | BioProject | Library | Tissues | Samples |
|---|---|---|---|---|---|
| Q. rubra multitissue | *Q. rubra* | PRJNA273270 | **Paired-end** | acorn, bud, leaf, twig, root, bark | 11 |
| Q. suber multitissue | *Q. suber* | PRJNA392919 | **Paired-end** | xylem, phellem, inner bark, leaf, pollen | 15 |
| Q. suber root | *Q. suber* | PRJNA690098 | **Single-end** | root tip, root segment | 8 |

## Directory Structure

```
02_gene_expression/
├── PE_paired_end/              ← For Q. rubra & Q. suber multitissue
│   ├── 0_downloadSRA_PE.sh
│   ├── 1_QualityCheck_PE.sh
│   ├── 2_Trimmomatic_PE.sh
│   └── 3_hisat2_mapping_PE.sh
├── SE_single_end/              ← For Q. suber root
│   ├── 0_downloadSRA_SE.sh
│   ├── 2_Trimmomatic_SE.sh
│   └── 3_hisat2_mapping_SE.sh
├── DEG_analysis.R              ← DESeq2 filtering for LBD genes
└── heatmap_visualization.R     ← Z-score heatmaps (pheatmap)
```

## PE vs SE Differences

The key differences between paired-end and single-end pipelines:

| Step | Paired-End (PE) | Single-End (SE) |
|---|---|---|
| SRA download | `fastq-dump --split-files` | `fastq-dump` (no split) |
| Trimmomatic | `PE` mode, 4 output files | `SE` mode, 1 output file |
| HISAT2 | `-1 R1 -2 R2` | `-U reads` |

## Setup

Replace `<YOUR_HPC_ID>` in all scripts with your HPC username before running.

## Acknowledgment

These scripts were adapted from the Functional Genomics course (BIOL7180) pipeline originally developed by **Dr. Rita Graze**, Auburn University.
