# CPDSeqer

This package is used to do CPD sequence data analysis.

# Prerequisites

Install Bowtie2 in Ubuntu with the following command: 

```
sudo apt-get install -y bowtie2
```

Install tabix in Ubuntu with the following command: 

```
sudo apt-get install –y tabix
```

Install Samtools in Ubuntu with the following command:

```
sudo apt-get install -y samtools
```

Install Python 3.7 in Ubuntu. Most factory versions of Ubuntu18.04 and later come with python pre-installed.  To check if Python is installed and the Python version, use the following command: 

```
python -version
```

If Python is not installed or version is lower than 3.7, use the following commands in sequential order to install Python 3.7:

```
sudo apt update
sudo apt install software-properties-common
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt update
sudo apt install python3.7
```

Install R packages in R

```
install.packages(c("knitr","rmarkdown"))
```

# Installation

Install CPD-seq protocol package in Ubuntu with the following command:

```
sudo pip3 install git+git://github.com/shengqh/cpdseqer.git
```

If you don't have pip installed, you need to install pip first.

```
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
sudo python get-pip.py
```

# Usage

## Demultiplex fastq file

If raw sequencing data in FASTQ format is multiplexed, it needs to be first de-multiplexed based on barcode sequence. Otherwise, skip to the next step. To de-multiplex a multiplexed FASTQ file, use the following command:

```
usage: cpdseqer demultiplex [-h] -i [INPUT] -o [OUTPUT] -b [BARCODE_FILE]

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT], --input [INPUT]
                        Input fastq file (gzipped supported)
  -o [OUTPUT], --output [OUTPUT]
                        Output folder
  -b [BARCODE_FILE], --barcode_file [BARCODE_FILE]
                        Tab-delimited file, first column is barcode, second column is sample name
```

for example:

```
cpdseqer demultiplex -i example.fastq.gz -o . -b barcode.txt
```

The [barcode.txt](https://cqsweb.app.vumc.org/download1/cpdseqer/data/barcode.txt) contains two columns indicate barcode and sample name (separated by tab).

|||
|---|---|
|ATCGCGAT|Control|
|GAACTGAT|UV|

You can download barcode and example.fastq.gz file:

```
wget https://cqsweb.app.vumc.org/download1/cpdseqer/data/example.fastq.gz
wget https://cqsweb.app.vumc.org/download1/cpdseqer/data/barcode.txt
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

## Count UV radiation induced DNA damage

This step completes raw data processing and generates important output files that will be required in multiple steps in the following workflow. 

```
usage: cpdseqer bam2dinucleotide [-h] -i [INPUT] -g [FASTA] [-q [MAPPING_QUALITY]] [-m MIN_COVERAGE] [-u] [-t] -o [OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT], --input [INPUT]
                        Input BAM file
  -g [FASTA], --fasta [FASTA]
                        Input genome fasta file
  -q [MAPPING_QUALITY], --mapping_quality [MAPPING_QUALITY]
                        Minimum mapping quality of read (default 20)
  -m MIN_COVERAGE, --min_coverage MIN_COVERAGE
                        The minimum coverage of dinucleotide for counting (default 1)
  -u, --unique_only     Use uniquely mapped read only
  -t, --test            Test the first 1000000 reads only
  -o [OUTPUT], --output [OUTPUT]
                        Output file prefix
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
wget https://cqsweb.app.vumc.org/download1/cpdseqer/data/UV.bam
wget https://cqsweb.app.vumc.org/download1/cpdseqer/data/UV.bam.bai
wget https://cqsweb.app.vumc.org/download1/cpdseqer/data/Control.bam
wget https://cqsweb.app.vumc.org/download1/cpdseqer/data/Control.bam.bai
```

## Quality control

QC based on dinucleotide count results can be performed using the following command:

```
usage: cpdseqer qc [-h] -i [INPUT] [-n [NAME]] [--count_type [COUNT_TYPE]] -o [OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT], --input [INPUT]
                        Input list file, first column is file name, second column is count file path, third column is dinucleotide file path
  -n [NAME], --name [NAME]
                        Input project name
  --count_type [COUNT_TYPE]
                        Input count type, rCnt/sCnt (read count/site count, default rCnt)
  -o [OUTPUT], --output [OUTPUT]
                        Output file prefix
```

for example:

```
cpdseqer qc -i dinucleotide.list -n qc -o qc_result
```

The [dinucleotide.list](https://cqsweb.app.vumc.org/download1/cpdseqer/data/dinucleotide.list) contains two columns indicate dinucleotide file and sample name (separated by tab).

|||
|---|---|
|Control.bed.bgz|Control|
|UV.dinucleotide.bed.bgz|UV|

You can download example files as following scripts.
```
wget https://cqsweb.app.vumc.org/download1/cpdseqer/data/UV.bed.bgz
wget https://cqsweb.app.vumc.org/download1/cpdseqer/data/UV.bed.bgz.tbi
wget https://cqsweb.app.vumc.org/download1/cpdseqer/data/UV.count
wget https://cqsweb.app.vumc.org/download1/cpdseqer/data/Control.bed.bgz
wget https://cqsweb.app.vumc.org/download1/cpdseqer/data/Control.bed.bgz.tbi
wget https://cqsweb.app.vumc.org/download1/cpdseqer/data/Control.count
wget https://cqsweb.app.vumc.org/download1/cpdseqer/data/dinucleotide.list
```

The [coordinates.list](https://cqsweb.app.vumc.org/download1/cpdseqer/coordinates.list) contains two columns indicate coordinate bed file and category name (separated by tab)

|||
|---|---|
|hg38_promoter.bed|Promoter|
|hg38_tf.bed|TFBinding|

You can download example files as following scripts. The hg38_promoter.bed contains three columns only. So, Promoter (from  coordinates.list definition) will be used as category name for all entries in the hg38_promoter.bed. The hg38_tf.bed contians four columns. The forth column in hg38_tf.bed indicates TF name which will be used as category name (--category_index 3) instead of TFBinding.

```
wget https://cqsweb.app.vumc.org/download1/cpdseqer/data/coordinates.list
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

# Running step by step

We provide example data and code for public download

```
wget -r --no-parent --reject "index.html*" https://cqsweb.app.vumc.org/download1/cpdseqer/
cd cqsweb.app.vumc.org/download1/cpdseqer/data
tar -xzvf GRCh38.p13.bowtie2.tar.gz
rm GRCh38.p13.bowtie2.tar.gz
cd ..
unzip test.zip
rm test.zip
```