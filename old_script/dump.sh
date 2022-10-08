cd fastq
mkdir dump
cell=$(awk 'NR!=1 {print $11}' ../ena.tsv | uniq)
echo $cell
echo $cell | while read -d " " id; do
	echo $id
	cat $(ls ${id}*.gz) >./dump/${id}.fastq.gz
done
