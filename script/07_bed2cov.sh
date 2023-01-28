#!/bin/bash
. $(which env_parallel.bash)
mkdir -p ../data/08_bed2cov/cov
# move the first column to the last column in all ../data/08_bed2cov/bed/*.bed files
function f_bed2cov() {
    file=$1
    cells=$(ls -1a ../data/04_post_align/*.bam | sed -E 's/\.\.\/data\/04_post_align\/[[:alnum:].]*_(.*)\.ucsc\.bam/\1/' | sort | uniq)
    # add .tmp to $file
    file_tmp=$file.tmp
    # move the first column to the last column
    awk -v FS='\t' -v OFS='\t' '{print $2,$3,$1}' $file >$file_tmp
    # get the 'chr' part in the third column, use regex, store it as the first column
    tmpc=$(awk -v FS='\t' -v OFS='\t' '{print $3}' $file_tmp | sed -E 's/.*_(chr[[:alnum:]]+)_.*/\1/')
    # insert tmpc as the first column
    paste -d '\t' <(echo "$tmpc") $file_tmp >$file.tmp2
    function bedtools_multicov() {
        cell=$1
        # get the *.bam files containing $cell, with relative path, use find, tab separated, sort
        bam=$(find ../data/04_post_align/ -name "*_$cell.ucsc.bam" | sort)
        # get the total count of reads in each bam file, use samtools view
        total_count=($(parallel samtools view -c ::: $bam))
        # divide each element in total_count by 1000000
        total_count=$(echo ${total_count[@]} | awk -v FS=' ' -v OFS=' ' '{for(i=1;i<=NF;i++){$i=$i/1000000};print}')
        # calculate cpm
        bedtools multicov -bams $bam -bed $file.tmp2 |
            awk -v FS='\t' -v OFS='\t' -v tc="$total_count" '{split(tc,tc2," ");print $1,$2,$3,$4,$5/tc2[1],$6/tc2[2],$7/tc2[3],$8/tc2[4]}' \
                >../data/08_bed2cov/cov/$(basename $file .JCEC.bed).$cell.cov
        # for sing_bam in $bam; do
        #     bedtools multicov -bams $sing_bam -b
        # done
    }
    export -f bedtools_multicov
    . $(which env_parallel.bash)
    env_parallel --env file bedtools_multicov ::: $cells
    rm $file_tmp && rm $file.tmp2
}
export -f f_bed2cov
parallel --progress f_bed2cov ::: ../data/08_bed2cov/bed/*.bed
