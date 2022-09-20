ls ./afterQC/*.gz | while read id; do bowtie2 -p 30 -x ./index/mm10 -U $id | samtools sort -O bam -o "./bam/${id:8: -14}.bam";done 
