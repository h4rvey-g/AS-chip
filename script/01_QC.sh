#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate fastp
# create output directory if not exist
mkdir -p ../data/02_fastq_afterQC
for i in ../data/01_fastq/*.fastq.gz; do
    echo "Processing $i"
    fastp -i $i \
        -o ../data/02_fastq_afterQC/$(basename $i .fastq.gz).fastq.gz \
        -h ../data/02_fastq_afterQC/report/$(basename $i .fastq.gz).html \
        -j ../data/02_fastq_afterQC/report/$(basename $i .fastq.gz).json \
        --phred64 \
        -w 8
done
# convert fastp json to MultiQC
multiqc ../data/02_fastq_afterQC/report/ -o ../data/02_fastq_afterQC/report/
