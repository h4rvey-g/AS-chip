conda init zsh
source ~/.zshrc
conda activate star
mkdir index bam
STAR --runThreadN 6 \
	--runMode genomeGenerate \
	--genomeDir ./index \
	--genomeFastaFiles ./Mus_musculus.GRCm39.dna.toplevel.fa \
	--sjdbGTFfile ./Mus_musculus.GRCm39.105.gtf \
	--sjdbOverhang 55
	IFS=" "
	echo $(ls ./fastq/dump/*.fastq.gz) | while read id; do
	#get cell name
	cell=$(echo $id|sed -E "s/.*dump\/(.*)\.fastq\.gz/\\1/")
	STAR --runThreadN 10 \
		--genomeDir ./index/ \
		--readFilesCommand zcat \
		--outSAMtype BAM SortedByCoordinate \
		--quantMode GeneCounts \
		--twopassMode Basic \
		--outFilterType BySJout \
		--outFilterMismatchNmax 999 \
		--outFilterMismatchNoverReadLmax 0.04 \
		--readFilesIn $id \
		--outFileNamePrefix ./bam/${cell};
	done
	cd ./bam
	# run rmats
	conda activate python
	ls *.bam|tr '[:space:]' ",">bam.txt
	~/.conda/envs/python/bin/python /home/guozhonghao/plugins/rmats-turbo/rmats.py --b1 ./bam.txt --bi ../index/ --gtf ../Mus_musculus.GRCm39.105.gtf --od ../output --tmp ../temp -t single --readLength 56 --libType fr-unstranded --nthread 10  --allow-clipping --statoff
