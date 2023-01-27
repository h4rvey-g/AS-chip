#!/bin/bash
for file in ../data/08_bed2cov/cov/*.cov; do
    # use sed to extract the site name from the file name
    site=$(echo $file | sed -E 's/.*\.(.*)\..*\.cov/\1/')
    if [ $site == "ExonStart" ]; then
        site="riExonStart"
    elif [ $site == "ExonEnd" ]; then
        site="riExonEnd"
    fi
    # get the 4-8 columns, tab separated
    awk -v FS='\t' -v OFS='\t' '{print $4,$5,$6,$7,$8}' $file |
        # remove site in the first column
        sed -E "s/_$site//" |
        #add column name, tab separated
        awk -v FS='\t' -v OFS='\t' -v st="$site" 'BEGIN{print "loci","H3K27Ac_"st,"H3K4me1_"st,"H3K4me2_"st,"H3K4me3_"st}{print $0}' >$file.tidy
done
cells=$(ls -1a ../data/04_post_align/*.bam | sed -E 's/.*_([[:alnum:]]*)\.ucsc\.bam/\1/' | sort | uniq)
as_type=(SE MXE A5SS A3SS RI)
for sing_as in ${as_type[@]}; do
    for sing_cell in $cells; do
        # get the file name
        files=$(find ../data/08_bed2cov/cov/ -name "$sing_as.*.$sing_cell.cov.tidy")
        #get the number of files
        num_files=$(echo $files | wc -w)
        # display num_files
        echo $num_files
        # paste the files, tab separated, remove columns with column name "loci"
        if [ $sing_as == "RI" ]; then
            paste -d '\t' $files | cut --complement -f6 >../data/08_bed2cov/$sing_as.$sing_cell.cov
        elif [ $num_files -eq 3 ]; then
            paste -d '\t' $files | cut --complement -f6,11 >../data/08_bed2cov/$sing_as.$sing_cell.cov
        elif [ $num_files -eq 4 ]; then
            paste -d '\t' $files | cut --complement -f6,11,16 >../data/08_bed2cov/$sing_as.$sing_cell.cov
        elif [ $num_files -eq 6 ]; then
            paste -d '\t' $files | cut --complement -f6,11,16,21,26 >../data/08_bed2cov/$sing_as.$sing_cell.cov
        fi
    done
done
