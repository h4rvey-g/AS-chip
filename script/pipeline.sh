# source ~/.zshrc
eval "$(conda shell.bash hook)"
conda activate /home/guozhonghao/.conda/envs/gcc/envs/chip-ap
chipap.py \
    --mode single \
    --ref mm10 \
    --peak broad \
    --genome /data0/work/guozhonghao/AS-chip/tools/ChIP-AP/chipap_installation/chipap_scripts/genomes \
    --output /data0/work/guozhonghao/AS-chip/data/output \
    --setname chip \
    --sample_table /data0/work/guozhonghao/AS-chip/data/sample.txt \
    # --custom_setting_table [path_to_setting_table_file].tsv \
    # --motif [path_to_known_motif_file] \
    --fcmerge --deltemp --run
