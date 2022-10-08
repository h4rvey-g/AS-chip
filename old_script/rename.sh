cd ./fastq/

awk 'NR!=1 {print}' ../ena.tsv | while read id; do
sra=$(echo $id | awk '{print $4}')
cell=$(echo $id | awk '{print $11}')

	if [[ -e $cell.fastq.gz ]]; then
		i=0
		while [[ -e $cell-$i.fastq.gz ]]; do
			let i++
		done
		cell=$cell-$i
	fi
	echo $sra
	mv $sra.fastq.gz $cell.fastq.gz
done
