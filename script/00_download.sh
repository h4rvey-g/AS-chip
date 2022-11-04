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
# batch rename files in the first column of sample.txt to the corresponding second column, begin at row 2
cat ../data/sample.txt | while read line; do
    # if $line begin with "sample", then skip
    if [[ $line =~ ^sample ]]; then
        continue
    fi
    old_name="../data/01_fastq/$(echo $line | cut -f 2)"
    new_name="../data/01_fastq/$(echo $line | cut -f 3)"
    mv $old_name $new_name
done
# traverse fastq.gz files in ../data/01_fastq/
for file in ../data/01_fastq/*.fastq.gz*~*; do
    old_name=$file
    new_name=$(echo $old_name | sed -r 's/(.*)\.fastq\.gz\.~([[:digit:]])~/\1_\2\.fastq\.gz/')
    mv $old_name $new_name
done
