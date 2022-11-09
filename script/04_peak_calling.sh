#!/bin/bash
# use macs2 to call peaks for bam in ../data/04_post_align/merged, output to ../data/05_peak_calling
mkdir -p ../data/05_peak_calling
find ../data/05_peak_calling/ -type f -delete
function call_peaks {
    macs2 callpeak -t $1 -f BAM -g mm -n $(basename $1 .ucsc_merged.bam) -B --broad --broad-cutoff 0.1 --extsize 200 --outdir ../data/05_peak_calling
}
export -f call_peaks
parallel --jobs 50 call_peaks ::: ../data/04_post_align/*.bam
