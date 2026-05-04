#!/bin/bash
#============================================================
#  Adapter Trimming & Quality Filtering (Trimmomatic) — PE
#
#  Removes sequencing adapters and low-quality bases from
#  paired-end reads. Produces paired and unpaired output.
#
#  Adapted from: Functional Genomics (BIOL7180) course scripts
#  Original author: Dr. Rita Graze, Auburn University
#  Modified by: Worlasi Dzamefe
#
#  HPC: Alabama Supercomputer Center (ASC)
#  Queue: medium | Cores: 6 | Time: 04:00:00 | Memory: 12gb
#============================================================

source /apps/profiles/modules_asax.sh.dyn
module load trimmomatic/0.39
module load fastqc/0.10.1

## Replace with your HPC username
MyID=<YOUR_HPC_ID>

WD=/scratch/${MyID}/RNAseq
DD=${WD}/RawData
CD=${WD}/CleanData
PCQ=${WD}/PostCleanQuality
adapters=AdaptersToTrim_All.fa

mkdir -p ${CD} ${PCQ}
cd ${DD}

## Build sample list from FASTQ filenames
ls | grep ".fastq" | cut -d "_" -f 1 | sort | uniq > list

## Copy adapter file (adjust path to your adapter FASTA)
cp /home/${MyID}/scripts/AdaptersToTrim_All.fa .

## Trim each sample (paired-end mode)
while read i
do
    java -jar /apps/x86-64/apps/spack_0.19.1/spack/opt/spack/linux-rocky8-zen3/gcc-11.3.0/trimmomatic-0.39-iu723m2xenra563gozbob6ansjnxmnfp/bin/trimmomatic-0.39.jar \
        PE -threads 6 -phred33 \
        ${DD}/${i}_1.fastq ${DD}/${i}_2.fastq \
        ${CD}/${i}_1_paired.fastq ${CD}/${i}_1_unpaired.fastq \
        ${CD}/${i}_2_paired.fastq ${CD}/${i}_2_unpaired.fastq \
        ILLUMINACLIP:${adapters}:2:30:10 \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:36

done < list

## Post-trim quality check
cd ${CD}
fastqc -t 6 -o ${PCQ} *_paired.fastq

cd ${WD}
tar -czf PostCleanQuality.tar.gz PostCleanQuality/
