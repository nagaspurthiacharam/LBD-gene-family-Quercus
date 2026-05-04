#!/bin/bash
#============================================================
#  Quality Assessment of Raw Paired-End Reads (FastQC)
#
#  Adapted from: Functional Genomics (BIOL7180) course scripts
#  Original author: Dr. Rita Graze, Auburn University
#  Modified by: Worlasi Dzamefe
#
#  HPC: Alabama Supercomputer Center (ASC)
#  Queue: medium | Cores: 1 | Time: 02:00:00 | Memory: 4gb
#============================================================

source /apps/profiles/modules_asax.sh.dyn
module load fastqc/0.10.1

## Replace with your HPC username
MyID=<YOUR_HPC_ID>

DD=/scratch/${MyID}/RNAseq/RawData
QC=/scratch/${MyID}/RNAseq/RawQuality

mkdir -p ${QC}
cd ${DD}

## Run FastQC on all raw FASTQ files
fastqc -t 4 -o ${QC} *.fastq

## Package results for transfer to local machine
cd /scratch/${MyID}/RNAseq
tar -czf RawQuality.tar.gz RawQuality/
