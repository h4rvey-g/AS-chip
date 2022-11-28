#!/bin/bash
# cat all .broadPeak files in ../data/05_peak_calling and use mergeBed to merge them
cat ../data/05_peak_calling/*.broadPeak | sort -k1,1 -k2,2n | mergeBed -i stdin >../data/07_merge_peaks/merged_peaks.broadPeak
# use bedtools multivoc to calculate reads abundance in certain ranges in bam files, and add title
bedtools multicov -bams ../data/04_post_align/*.bam -bed ../data/07_merge_peaks/merged_peaks.broadPeak >../data/07_merge_peaks/merged_peaks.broadPeak.cov
bam_names=$(ls ../data/04_post_align/*.bam | sed 's/.*\///' | sed 's/.bam//' | tr "\n" "\t")
echo -e "chr\tstart\tend\t$bam_names" >../data/07_merge_peaks/merged_peaks.broadPeak.cov.tmp
cat ../data/07_merge_peaks/merged_peaks.broadPeak.cov >>../data/07_merge_peaks/merged_peaks.broadPeak.cov.tmp
mv ../data/07_merge_peaks/merged_peaks.broadPeak.cov.tmp ../data/07_merge_peaks/merged_peaks.broadPeak.cov
