# CPDSeqer

This package is used to do CPD sequence data analysis.

# Prerequisites

Install tabix

```
apt-get install -y tabix
```

# Installation

Install python main package

```
pip install git+git://github.com/shengqh/cpdseqer.git
```

# Usage

## Demultiplex fastq file

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

Example of [barcode.txt](https://github.com/shengqh/cpdseqer/raw/master/data/barcode.txt), columns are separated by tab
```
ATCGCGAT        Control
GAACTGAT        UV
```

## Align reads to genome using bowtie2

```
bowtie2 -p 8  -x hg38/bowtie2_index_2.3.5.1/GRCh38.p12.genome -U Control.fastq.gz \
  --sam-RG ID:Control --sam-RG LB:Control --sam-RG SM:Control --sam-RG PL:ILLUMINA \
  -S Control.sam 2> Control.log
  
samtools view -Shu -F 256 Control.sam | samtools sort -o Control.bam -T Control -

samtools index Control.bam

rm Control.sam
```

You can download hg38 bowtie2 index files:

```
mkdir hg38
cd hg38
wget https://cqsweb.app.vumc.org/download1/cpdseqer/hg38_bowtie2.tar.gz
tar -xzvf hg38_bowtie2.tar.gz
cd ..
```

## Extract dinucleotide from bam file

```
usage: cpdseqer bam2dinucleotide [-h] -i [INPUT] -g [GENOME_SEQ_FILE]
                                 [-q [MAPPING_QUALITY]] -o [OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT], --input [INPUT]
                        Input BAM file
  -g [GENOME_SEQ_FILE], --genome_seq_file [GENOME_SEQ_FILE]
                        Input genome seq file
  -q [MAPPING_QUALITY], --mapping_quality [MAPPING_QUALITY]
                        Minimum mapping quality of read (default 20)
  -o [OUTPUT], --output [OUTPUT]
                        Output file name
```

for example:

```
cpdseqer bam2dinucleotide \
  -i Control.bam \
  -g hg38/bowtie2_index_2.3.4.1/GRCh38.p12.genome.fa \
  -o Control.dinucleotide.bed.bgz
```

You can download bam files:

```
wget https://cqsweb.app.vumc.org/download1/cpdseqer/UV.bam
wget https://cqsweb.app.vumc.org/download1/cpdseqer/UV.bam.bai
wget https://cqsweb.app.vumc.org/download1/cpdseqer/Control.bam
wget https://cqsweb.app.vumc.org/download1/cpdseqer/Control.bam.bai
```

## Get statistic based on bed file

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
cpdseqer statistic --category_index 3 \
  -i cpd__fileList1.list \
  -c cpd__fileList2.list \
  -o cpd.dinucleotide.statistic.tsv
```

cpd__fileList1.list (separated by tab)

```
Control.dinucleotide.bed.bgz      Control
UV.dinucleotide.bed.bgz           UV
```

cpd__fileList2.list (separated by tab)

```
hg38_promoter.bed      Promoter
hg38_tf.bed            TFBinding
```

You can download example files as following scripts. The hg38_promoter.bed contains three columns only. So, Promoter (from  cpd__fileList2.list definition) will be used as category name for all entries in the hg38_promoter.bed. The hg38_tf.bed contians four columns. The forth column in hg38_tf.bed indicates TF name which will be used as category name (--category_index 3) instead of TFBinding.

```
wget https://cqsweb.app.vumc.org/download1/cpdseqer/UV.dinucleotide.bed.bgz
wget https://cqsweb.app.vumc.org/download1/cpdseqer/UV.dinucleotide.bed.bgz.tbi
wget https://cqsweb.app.vumc.org/download1/cpdseqer/Control.dinucleotide.bed.bgz
wget https://cqsweb.app.vumc.org/download1/cpdseqer/Control.dinucleotide.bed.bgz.tbi
wget https://cqsweb.app.vumc.org/download1/cpdseqer/cpd__fileList1.list
wget https://cqsweb.app.vumc.org/download1/cpdseqer/cpd__fileList2.list
wget https://github.com/shengqh/cpdseqer/raw/master/data/hg38_promoter.bed
wget https://github.com/shengqh/cpdseqer/raw/master/data/hg38_tf.bed


```

# Running cpdseqer using singularity

We also build docker container for cpdseqer which can be used by singularity.

## Running directly

```
singularity exec -e docker://shengqh/cpdseqer bowtie2 -h
singularity exec -e docker://shengqh/cpdseqer cpdseqer -h
```

## Convert docker image to singularity image first

```
singularity build cpdseqer.simg docker://shengqh/cpdseqer
singularity exec -e cpdseqer.simg bowtie2 -h
singularity exec -e cpdseqer.simg cpdseqer -h
```
