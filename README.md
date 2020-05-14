# CPDSeqer

This package is used to do CPD sequence data analysis.

# Prerequisites

Install tabix, samtools and bowtie2

```
sudo apt-get install -y tabix samtools bowtie2 pandoc
```

Install R packages

```
install.packages(c("knitr","rmarkdown"))
```

# Installation

Install cpdseqer package from github

```
sudo pip install git+git://github.com/shengqh/cpdseqer.git
```

If you don't have pip installed, you need to install pip first.

```
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
sudo python get-pip.py
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
cpdseqer demultiplex -i example.fastq.gz -o . -b barcode.txt
```

The [barcode.txt](https://cqsweb.app.vumc.org/download1/cpdseqer/barcode.txt) contains two columns indicate barcode and sample name (separated by tab).

|||
|---|---|
|ATCGCGAT|Control|
|GAACTGAT|UV|

You can download barcode and example.fastq.gz file:

```
wget https://cqsweb.app.vumc.org/download1/cpdseqer/example.fastq.gz
wget https://cqsweb.app.vumc.org/download1/cpdseqer/barcode.txt
```

## Align reads to genome using bowtie2

```
bowtie2 -p 8  -x hg38/bowtie2_index_2.3.5.1/GRCh38.p12.genome -U Control.fastq.gz -S Control.sam 2> Control.log
  
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
                                 [-q [MAPPING_QUALITY]] [-u] -o [OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT], --input [INPUT]
                        Input BAM file
  -g [GENOME_SEQ_FILE], --genome_seq_file [GENOME_SEQ_FILE]
                        Input genome seq file
  -q [MAPPING_QUALITY], --mapping_quality [MAPPING_QUALITY]
                        Minimum mapping quality of read (default 20)
  -u, --unique_only     Use uniquely mapped read only
  -o [OUTPUT], --output [OUTPUT]
                        Output file name

```

for example:

```
cpdseqer bam2dinucleotide \
  -i Control.bam \
  -g hg38/bowtie2_index_2.3.4.1/GRCh38.p12.genome.fa \
  -o Control
```

You can download our example bam files:

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
  -i dinucleotide.list \
  -c coordinates.list \
  -o cpd.dinucleotide.statistic.tsv
```

The [dinucleotide.list](https://cqsweb.app.vumc.org/download1/cpdseqer/dinucleotide.list) contains two columns indicate dinucleotide file and sample name (separated by tab).

|||
|---|---|
|Control.bed.bgz|Control|
|UV.dinucleotide.bed.bgz|UV|

The [coordinates.list](https://cqsweb.app.vumc.org/download1/cpdseqer/coordinates.list) contains two columns indicate coordinate bed file and category name (separated by tab)

|||
|---|---|
|hg38_promoter.bed|Promoter|
|hg38_tf.bed|TFBinding|

You can download example files as following scripts. The hg38_promoter.bed contains three columns only. So, Promoter (from  coordinates.list definition) will be used as category name for all entries in the hg38_promoter.bed. The hg38_tf.bed contians four columns. The forth column in hg38_tf.bed indicates TF name which will be used as category name (--category_index 3) instead of TFBinding.

```
wget https://cqsweb.app.vumc.org/download1/cpdseqer/UV.bed.bgz
wget https://cqsweb.app.vumc.org/download1/cpdseqer/UV.bed.bgz.tbi
wget https://cqsweb.app.vumc.org/download1/cpdseqer/UV.count
wget https://cqsweb.app.vumc.org/download1/cpdseqer/Control.bed.bgz
wget https://cqsweb.app.vumc.org/download1/cpdseqer/Control.bed.bgz.tbi
wget https://cqsweb.app.vumc.org/download1/cpdseqer/Control.count
wget https://cqsweb.app.vumc.org/download1/cpdseqer/dinucleotide.list
wget https://cqsweb.app.vumc.org/download1/cpdseqer/coordinates.list
wget https://github.com/shengqh/cpdseqer/raw/master/data/hg38_promoter.bed
wget https://github.com/shengqh/cpdseqer/raw/master/data/hg38_tf.bed

```

## Draw position figure

```
usage: cpdseqer fig_position [-h] -i [INPUT] -c [COORDINATE_FILE] [-b [BACKGROUND_FILE]] [-s] [--add_chr] [-t] -o [OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT], --input [INPUT]
                        Input dinucleotide list file, first column is file location, second column is file name
  -c [COORDINATE_FILE], --coordinate_file [COORDINATE_FILE]
                        Input coordinate bed file (can use short name hg38/hg19 as default nucleosome file)
  -b [BACKGROUND_FILE], --background_file [BACKGROUND_FILE]
                        Background dinucleotide file
  -s, --space           Use space rather than tab in coordinate files
  --add_chr             Add chr in chromosome name in coordinate file
  -t, --test            Test the first 10000 coordinates only
  -o [OUTPUT], --output [OUTPUT]
                        Output file name
```

for example, we will calculate the dinucleotide position in nucleosome:

```
cpdseqer fig_position -s -b hg38_background.bed.bgz -i dinucleotide.list -c hg38 -o cpd_position.txt
```

Here, you can input absolute coordinate file, or hg38/hg19. hg38 and hg19 indicates the nucleosome coordinate files which can be downloaded by:

```
wget https://github.com/shengqh/cpdseqer/raw/master/cpdseqer/data/nucleosome_hg19_interval.zip
wget https://github.com/shengqh/cpdseqer/raw/master/cpdseqer/data/nucleosome_hg38_interval.zip

```

You can download hg38 and yeast background file from:

```
wget https://cqsweb.app.vumc.org/download1/cpdseqer/XXXXX.bed.bgz
wget https://cqsweb.app.vumc.org/download1/cpdseqer/XXXXX.bed.bgz.tbi
wget https://cqsweb.app.vumc.org/download1/cpdseqer/XXXXX.count
wget https://cqsweb.app.vumc.org/download1/cpdseqer/XXXXX.bed.bgz
wget https://cqsweb.app.vumc.org/download1/cpdseqer/XXXXX.bed.bgz.tbi
wget https://cqsweb.app.vumc.org/download1/cpdseqer/XXXXX.count
```

## Draw genome/chromosome figure

```
usage: cpdseqer fig_genome [-h] -i [INPUT] [-b [BLOCK]] [-d [DB]] -o [OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT], --input [INPUT]
                        Input configuration file
  -b [BLOCK], --block [BLOCK]
                        Block size for summerize dinucleotide count
  -d [DB], --db [DB]    Input database version, hg38 or hg19, default is hg38
  -o [OUTPUT], --output [OUTPUT]
                        Output file name
```

for example:

```
cpdseqer fig_genome -i fig_genome.input.txt -d hg38 -o cpd.report
```

[fig_genome.input.txt](https://github.com/shengqh/cpdseqer/raw/master/cpdseqer/data/fig_genome.input.txt) contains three columns indicate name, file and group (1=case and 0=control) which are seprated by tab.

|Name|File|Case|
|---|---|---|
|Control|Control.bed.bgz|0|
|UV|UV.bed.bgz|1|

You can download it from:

```
wget https://github.com/shengqh/cpdseqer/raw/master/cpdseqer/data/fig_genome.input.txt
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
