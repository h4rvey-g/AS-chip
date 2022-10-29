#!/bin/bash
# if directory index does not exist, run bwa build index
if [ ! -d ../data/genome/index ]; then
    mkdir -p ../data/genome/index
    bwa index -p ../data/genome/index/GRCm39 ../data/genome/Mus_musculus.GRCm39.dna.toplevel.fa
fi
# Use GNU parallel to run bwa mem
mkdir -p ../data/03_align
function BWA_align() {
    echo "Processing $1"
    # First unzip fastq.gz to fq, then use BWA mem to align, and convert to bam
    gunzip -c $1 |
        bwa mem -t 2 -M ../data/genome/index/GRCm39 - |
        samtools view -@ 2 -bS - \
            >../data/03_align/$(basename $1 .fastq.gz).bam
}
export -f BWA_align
parallel --jobs 50 BWA_align ::: ../data/02_fastq_afterQC/*.fastq.gz
