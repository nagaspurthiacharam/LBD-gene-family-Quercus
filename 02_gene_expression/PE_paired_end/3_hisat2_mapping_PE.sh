#!/bin/bash
#============================================================
#  Read Mapping with HISAT2 & Transcript Quantification — PE
#
#  Maps paired-end reads to reference genome, then assembles
#  and quantifies transcripts with StringTie.
#
#  Adapted from: Functional Genomics (BIOL7180) course scripts
#  Original author: Dr. Rita Graze, Auburn University
#  Modified by: Worlasi Dzamefe
#
#  HPC: Alabama Supercomputer Center (ASC)
#  Queue: medium | Cores: 6 | Time: 06:00:00 | Memory: 16gb
#============================================================

source /apps/profiles/modules_asax.sh.dyn
module load hisat2/2.2.0
module load stringtie/2.2.1
module load samtools
module load gcc
module load python/3.10.8-zimemtc

## Replace with your HPC username
MyID=<YOUR_HPC_ID>

WD=/scratch/${MyID}/RNAseq
CD=${WD}/CleanData
RONE=_1_paired.fastq
RTWO=_2_paired.fastq
RONE_UN=_1_unpaired.fastq
RTWO_UN=_2_unpaired.fastq

## Reference genome and annotation
## Adjust these paths to your genome files
REF=/home/${MyID}/genomes/genome
GFF=/home/${MyID}/genomes/genome_annotation.gff3

## Build genome index (run once, then comment out)
# hisat2-build ${REF}.fa ${REF}

## Mapping output directory
MAPDIR=${WD}/Map_HiSat2
mkdir -p ${MAPDIR}

cd ${CD}
ls | grep "_1_paired.fastq" | cut -d "_" -f 1 | sort | uniq > list

while read i
do
    ## Map paired-end reads with HISAT2
    ## --dta: report alignments tailored for transcript assemblers
    hisat2 -p 6 --dta \
        -x ${REF} \
        -1 ${CD}/${i}${RONE} \
        -2 ${CD}/${i}${RTWO} \
        -S ${MAPDIR}/${i}.sam

    ## Convert SAM to sorted BAM
    samtools sort -@ 6 -o ${MAPDIR}/${i}_sorted.bam ${MAPDIR}/${i}.sam
    samtools index ${MAPDIR}/${i}_sorted.bam

    ## Remove SAM to save space
    rm ${MAPDIR}/${i}.sam

    ## Assemble and quantify transcripts with StringTie
    stringtie -p 6 -e -G ${GFF} \
        -o ${MAPDIR}/${i}_stringtie.gtf \
        -A ${MAPDIR}/${i}_gene_abundance.tab \
        ${MAPDIR}/${i}_sorted.bam

done < list

## Generate count matrix for DESeq2
## Create sample list file: sample_id  path_to_gtf
cd ${MAPDIR}
for f in *_stringtie.gtf; do
    sample=$(echo $f | sed 's/_stringtie.gtf//')
    echo "${sample} ${MAPDIR}/${f}"
done > sample_list.txt

## Extract count matrix using prepDE.py (provided with StringTie)
python3 prepDE.py -i sample_list.txt -g gene_count_matrix.csv -t transcript_count_matrix.csv
