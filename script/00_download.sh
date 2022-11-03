#!/bin/bash
# use fasterq-dump to download fastq files in accession.txt
# if files exist in ../data/01_fastq, let user choose to delete theme or not
if [ -d ../data/01_fastq ]; then
    read -p "Files exist in ../data/01_fastq, do you want to delete them? [y/n] " -n 1 -r REPLY
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        rm -rf ../data/01_fastq/*
    fi
else
    mkdir -p ../data/01_fastq
fi
function getFastq() {
    echo "Processing $1"
    # use fasterq-dump to download fastq files
    fasterq-dump \
        --split-files \
        --skip-technical \
        --threads 2 \
        --temp ../data/01_fastq \
        $1 -O ../data/01_fastq
}
export -f getFastq
parallel getFastq ::: $(cat ../data/accession.txt)
# use parallel to gzip compress fastq files, and delete the original files
parallel gzip {} ::: ../data/01_fastq/*.fastq
# use fastqc to check the quality of fastq files, store the results in ../data/01_fastq/report_fastqc, and use multiqc to generate a report
mkdir -p ../data/01_fastq/report_fastqc
fastqc -o ../data/01_fastq/report_fastqc ../data/01_fastq/*.fastq.gz
/data0/apps/anaconda3/bin/multiqc ../data/01_fastq/report_fastqc/ -o ../data/01_fastq/report_fastqc/
