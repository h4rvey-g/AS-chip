#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate fastp
# create output directory if not exist
mkdir -p ../data/02_fastq_afterQC && mkdir -p ../data/02_fastq_afterQC/report
function getQC() {
    echo "Processing $1"
    fastp -i $1 \
        -o ../data/02_fastq_afterQC/$(basename $1 .fastq.gz).fastq.gz \
        -h ../data/02_fastq_afterQC/report/$(basename $1 .fastq.gz).html \
        -j ../data/02_fastq_afterQC/report/$(basename $1 .fastq.gz).json \
        --phred64 \
        -f 6 \
        -Q
}
export -f getQC
parallel getQC ::: ../data/01_fastq/*.fastq.gz
# convert fastp json to MultiQC
/data0/apps/anaconda3/bin/multiqc ../data/02_fastq_afterQC/report/ \
    -o ../data/02_fastq_afterQC/report/ \
    -f \
    -c ../data/02_fastq_afterQC/report/multiqc_config.yaml
