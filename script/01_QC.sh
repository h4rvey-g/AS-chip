#!/bin/bash
# eval "$(conda shell.bash hook)"
# conda activate fastp
# create output directory if not exist
mkdir -p ../data/02_fastq_afterQC && mkdir -p ../data/02_fastq_afterQC/report
function getQC() {
    echo "Processing $1"
    # use trim_galore to remove adapter, and use trimmomatic to cut first 10bp
    trim_galore \
        --phred64 \
        --gzip \
        --illumina \
        --clip_R1 10 \
        --length 20 \
        --stringency 3 \
        --fastqc --fastqc_args "-o ../data/02_fastq_afterQC/report" \
        -o ../data/02_fastq_afterQC $1
}
export -f getQC
parallel getQC ::: ../data/01_fastq/*.fastq.gz
# convert fastp json to MultiQC
/data0/apps/anaconda3/bin/multiqc ../data/02_fastq_afterQC/report/ \
    -o ../data/02_fastq_afterQC/report/ \
    -f