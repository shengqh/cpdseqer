
cd /scratch/cqs/shengq2/guoyan/test/uv_comp_genome_region

#python /home/shengq2/program/cpdseqer/debug.py uv_comp_genome_region -i dinucleotide.list -d hg38.1000000.fa -c nucleosome_hg38_1000.txt -o uv_comp_genome -s
python /home/shengq2/program/cpdseqer/debug.py uv_comp_genome_region -i dinucleotide.list -d hg38.1000000.fa -c nucleosome_hg38_interval.txt -o uv_comp_genome -s
