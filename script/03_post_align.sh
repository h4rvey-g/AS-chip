#!/bin/bash
# Use samtools to sort and index bam files in 03_align, and store in 04_post_align
mkdir -p ../data/04_post_align
function samtools_sort_index() {
    echo "Processing $1"
    samtools sort -@ 2 -o ../data/04_post_align/$(basename $1 .bam)_sorted.bam $1
    samtools index -@ 2 ../data/04_post_align/$(basename $1 .bam)_sorted.bam
}
export -f samtools_sort_index
parallel --jobs 50 samtools_sort_index ::: ../data/03_align/*.bam
# Use fastqc to check the quality of the bam files, store results in 04_post_align/report, and use multiqc to summarize the results
mkdir -p ../data/04_post_align/report
fastqc -t 64 -o ../data/04_post_align/report ../data/04_post_align/*.bam
multiqc ../data/04_post_align/report -o ../data/04_post_align/report
# Use a for loop to traverse sample names in sample_list.txt, and use samtools to merge bam files for each sample
mkdir -p ../data/04_post_align/merged
for sample in $(cat ../data/sample_list.txt); do
    echo "Processing $sample"
    samtools merge -@ 50 ../data/04_post_align/merged/${sample}_merged.bam ../data/04_post_align/${sample}*.bam
done