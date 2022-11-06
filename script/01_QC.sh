#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate fastx
# create output directory if not exist
mkdir -p ../data/02_fastq_afterQC && mkdir -p ../data/02_fastq_afterQC/report
find ../data/02_fastq_afterQC/ -type f -delete
function getQC() {
    echo "Processing $1"
    # first unzip, then use fastx_trimmer to trim the first 10bp
    gunzip -c $1 | fastx_trimmer -f 10 -z -o ../data/02_fastq_afterQC/$(basename $1)
    # perform fastqc
    fastqc ../data/02_fastq_afterQC/$(basename $1) -o ../data/02_fastq_afterQC/report
}
export -f getQC
parallel getQC ::: ../data/01_fastq/*.fastq.gz
# convert fastp json to MultiQC
/data0/apps/anaconda3/bin/multiqc ../data/02_fastq_afterQC/report/ \
    -o ../data/02_fastq_afterQC/report/ \
    -f
