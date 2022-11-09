#!/bin/bash
# Use samtools to sort and index bam files in 03_align, and store in 04_post_align
mkdir -p ../data/04_post_align
mkdir -p ../data/04_post_align/merged
find ../data/04_post_align/ -type f -delete
function sort_index {
    echo "Processing $1"
    samtools merge -@ 4 -f ../data/04_post_align/$1.bam ../data/03_align/$1*.bam
}
export -f sort_index
cat ../data/sample_list.txt | parallel --jobs 50 sort_index
function samtools_sort_index() {
    echo "Processing $1"
    samtools sort -@ 2 -o $1 $1
    samtools index -@ 2 $1
}
export -f samtools_sort_index
parallel --jobs 50 samtools_sort_index ::: ../data/04_post_align/*.bam
# Use fastqc to check the quality of the bam files, store results in 04_post_align/report, and use multiqc to summarize the results
mkdir -p ../data/04_post_align/report
fastqc -t 64 -o ../data/04_post_align/report ../data/04_post_align/*.bam
multiqc ../data/04_post_align/report -o ../data/04_post_align/report -f
