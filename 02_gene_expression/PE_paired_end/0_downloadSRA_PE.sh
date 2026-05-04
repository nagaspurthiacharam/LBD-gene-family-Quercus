#!/bin/bash
#============================================================
#  Download Paired-End RNA-seq Data from NCBI SRA
#
#  Datasets:
#    Q. rubra multitissue — PRJNA273270 (11 samples, 6 tissues)
#    Q. suber multitissue — PRJNA392919 (15 samples, 5 tissues)
#
#  Adapted from: Functional Genomics (BIOL7180) course scripts
#  Original author: Dr. Rita Graze, Auburn University
#  Modified by: Worlasi Dzamefe
#
#  HPC: Alabama Supercomputer Center (ASC)
#  Queue: medium | Cores: 1 | Time: 04:00:00 | Memory: 4gb
#============================================================

########## Load Modules
source /apps/profiles/modules_asax.sh.dyn
module load sra

########## Define variables
## Replace MyID with your HPC username
MyID=<YOUR_HPC_ID>               ## e.g., MyID=aubxyz001

## Working and data directories (adjust paths as needed)
WD=/scratch/${MyID}/RNAseq
DD=${WD}/RawData

mkdir -p ${DD}
cd ${DD}

## Force SRA config initialization
vdb-config --interactive > /dev/null 2>&1 <<EOF
q
EOF

##########  Q. rubra — PRJNA273270 (Paired-End)
## 11 samples across 6 tissues (2 genotypes: Father, Mother)
## Tissues: acorn, bud, immature leaf, mature leaf, immature twig, root
fastq-dump -F --split-files SRR1778917
fastq-dump -F --split-files SRR1778919
fastq-dump -F --split-files SRR1778921
fastq-dump -F --split-files SRR1778922
fastq-dump -F --split-files SRR1778923
fastq-dump -F --split-files SRR1778925
fastq-dump -F --split-files SRR1778927
fastq-dump -F --split-files SRR1778929
fastq-dump -F --split-files SRR1778930
fastq-dump -F --split-files SRR1778931
fastq-dump -F --split-files SRR1778932

##########  Q. suber multitissue — PRJNA392919 (Paired-End)
## 15 samples across 5 tissues: xylem, phellem, inner bark, leaf, pollen
## Add Q. suber SRA accessions below (from NCBI BioProject PRJNA392919)
# fastq-dump -F --split-files SRR_ACCESSION1
# fastq-dump -F --split-files SRR_ACCESSION2
# ... (add all accessions from the BioProject)

## --split-files: splits paired-end reads into _1.fastq and _2.fastq
## -F: defline contains only original sequence name
