# cpdseqer

This package is used to do CPD sequence data analysis.

# installation

Install python main package

```
pip install git+git://github.com/shengqh/cpdseqer.git
```

Install tabix

```
apt-get install -y tabix
```

# usage

### Demultiplex fastq file

```
usage: cpdseqer demultiplex [-h] -i [INPUT] -o [OUTPUT] -b [BARCODEFILE]

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT], --input [INPUT]
                        Input fastq gzipped file
  -o [OUTPUT], --output [OUTPUT]
                        Output folder
  -b [BARCODEFILE], --barcodeFile [BARCODEFILE]
                        Tab-delimited file, first column is barcode, second
                        column is sample name
```
for example:
```
cpdseqer demultiplex -i sample.fastq.gz -o . -b barcode.txt
```

barcode.txt (separated by tab)
```
ATCGCGAT        Control
GAACTGAT        UV
```

### Align reads to genome using bowtie2

```
bowtie2 -p 8  -x hg38/bowtie2_index_2.3.4.1/GRCh38.p12.genome -U Control.fastq.gz --sam-R
G ID:Control --sam-RG LB:Control --sam-RG SM:Control --sam-RG PL:ILLUMINA -S Control.sam 2> Control.log
samtools view -Shu -F 256 Control.sam | samtools sort -o Control.bam -T Control -
samtools index Control.bam
rm Control.sam
```

### Extract dinucleotide from bam file

```
usage: cpdseqer bam2dinucleotide [-h] -i [INPUT] -g [GENOME_SEQ_FILE] -o
                                 [OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT], --input [INPUT]
                        Input BAM file
  -g [GENOME_SEQ_FILE], --genome_seq_file [GENOME_SEQ_FILE]
                        Input genome seq file
  -o [OUTPUT], --output [OUTPUT]
                        Output file name
```

for example:

```
cpdseqer bam2dinucleotide -i Control.bam -g hg38/bowtie2_index_2.3.4.1/GRCh38.p12.genome.fa -o Control.dinucleotide.bed.bgz
```

### Get statistic based on bed file

```
usage: cpdseqer statistic [-h] -i [INPUT] -c [COORDINATE_LIST_FILE] [-s]
                          [--category_index CATEGORY_INDEX] [--add_chr] -o
                          [OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT], --input [INPUT]
                        Input dinucleotide list file, first column is file
                        location, second column is file name
  -c [COORDINATE_LIST_FILE], --coordinate_list_file [COORDINATE_LIST_FILE]
                        Input coordinate list file, first column is file
                        location, second column is file name
  -s, --space           Use space rather than tab in coordinate files
  --category_index CATEGORY_INDEX
                        Zero-based category column index in coordinate file
  --add_chr             Add chr in chromosome name in coordinate file
  -o [OUTPUT], --output [OUTPUT]
                        Output file name

```

for example:

```
cpdseqer statistic --category_index 6 -i cpd__fileList1.list -c cpd__fileList2.list -o cpd.dinucleotide.statistic.tsv
```

cpd__fileList1.list (separated by tab)

```
/scratch/dinucleotide/Control.dinucleotide.bed.bgz      Control
/scratch/dinucleotide/UV.dinucleotide.bed.bgz           UV
```

cpd__fileList2.list (separated by tab)

```
/scratch/hg38_promoter.bed      Promoter
/scratch/hg38_tf.bed            TFBinding
```


# running cpdseqer using singularity

We also build docker container for cpdseqer which can be used by singularity.

## running directly

```
singularity exec -e docker://shengqh/cpdseqer cpdseqer -h
```

## convert docker image to singularity image first

```
singularity build cpdseqer.simg docker://shengqh/cpdseqer
singularity exec -e cpdseqer.simg cpdseqer -h
```
