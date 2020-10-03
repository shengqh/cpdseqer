# CPDSeqer

This package is used to do CPD sequence data analysis.

# Table of Contents
1. [Prerequisites](#Prerequisites)
2. [Installation](#Installation)
3. [Usage](#Usage)
    1. [(Optional) Demultiplex fastq file](#1-optional-demultiplex-fastq-file)
    2. [General QC](#2-general-qc)
    3. [Generate index files for reference genome (FASTA) files](#3-generate-index-files-for-reference-genome-fasta-files)
    4. [Align reads to genome using bowtie2](#4-align-reads-to-genome-using-bowtie2)
    5. [(Optional) Correct GC content bias](#5-optional-correct-gc-content-bias)
    6. [Count UV radiation induced DNA damage](#6-count-uv-radiation-induced-dna-damage)
    7. [(Optional) Generate text files to inform estimated lesions along chromosomes](#7-optional-generate-text-files-to-inform-estimated-lesions-along-chromosomes)
    8. [(Optional) Subtract out short tandem repeat regions or narrow down to the interested gnomic regions](#8-optional-subtract-out-short-tandem-repeat-regions-or-narrow-down-to-the-interested-gnomic-regions)
    9. [(Optional) Estimate sample-wise normalization factors](#9-optional-estimate-sample-wise-normalization-factors)
    10. [Quality control](#10-quality-control)
    11. [Generate a genome-wide UV damage distribution map](#11-generate-a-genome-wide-uv-damage-distribution-map)
    12. [Draw dinucleotide pileup figure in a specific genomic region type](#12-draw-dinucleotide-pileup-figure-in-a-specific-genomic-region-type)
    13. [Compare UV radiation damage of sample(s) against the reference genome background](#13-compare-uv-radiation-damage-of-samples-against-the-reference-genome-background)
    14. [Compare UV damage of sample(s) against reference genome background within a specific region type](#14-compare-uv-damage-of-samples-against-reference-genome-background-within-a--specific-region-type)
    15. [Compare UV damage between two regions for one or multiple samples](#15-compare-uv-damage-between-two-regions-for-one-or-multiple-samples)
    16. [Compare genome-wide UV damage between two groups of samples](#16-compare-genome-wide-uv-damage-between-two-groups-of-samples)
    17. [Compare UV damage between two groups of samples within a specific region type](#17-compare-uv-damage-between-two-groups-of-samples-within-a-specific-region-type)
4. [Running cpdseqer using singularity](#running-cpdseqer-using-singularity)

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

Install bedtools in Ubuntu with the following command:

```
sudo apt-get install -y bedtools
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
install.packages(c("knitr","rmarkdown", "data.table", "R.utils", "ggplot2", "reshape2"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")
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

## 1. (Optional) Demultiplex fastq file

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

The [barcode.txt](https://cqsweb.app.vumc.org/Data/cpdseqer/data/barcode.txt) contains two columns indicate barcode and sample name (separated by tab).

|||
|---|---|
|ATCGCGAT|Control|
|GAACTGAT|UV|

You can download barcode and example.fastq.gz file:

```
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/example.fastq.gz
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/barcode.txt
```

An alternative method to perform demultiplex is using [Je](https://gbcs.embl.de/portal/tiki-index.php?page=Je).

```
je demultiplex F1=[INPUT] BF=[JE_BARCODEFILE] O=[OUTPUT]
```

The [JE_BARCODEFILE](https://cqsweb.app.vumc.org/Data/cpdseqer/data/barcode_je.txt) contains three columns indicate sample name, barcode and sample filename (separated by tab).

||||
|---|---|---|
|Control|ATCGCGAT|Control.fastq.gz|
|UV|GAACTGAT|UV.fastq.gz|

## 2. General QC
Without going into detail, we recommend previously established methods such as FASTQC and QC3 for this general QC step.

## 3. Generate index files for reference genome (FASTA) files
```
bowtie2-build [–-threads [THREADS]][INPUT] [OUTPUT]
```
## 4. Align reads to genome using bowtie2
(A)	If sequencing data is single-end, use the following command:

```
bowtie2 -p [THREADS] -x [INDEX] -U [FASTQ] -S [SAM] | samtools sort –o [OUTPUT] –T [TEMP_PREFIX] [-@ [THREADS]] [–m [MAX_MEMORY]]
```
(B)	If sequencing data is pair-end, use the following command:
```
bowtie2 -p [THREADS] -x [INDEX] -1 [FASTQ1] -2 [FASTQ2] –S [SAM] | samtools sort –o [OUTPUT] –T [TEMP_PREFIX] [-@ [THREADS]] [–m [MAX_MEMORY]]
```
for example,
```
bowtie2 -p 8  -x hg38/bowtie2_index_2.3.5.1/GRCh38.p12.genome -U Control.fastq.gz | samtools sort -@ 8 -m 4G -o Control.bam -T Control -
samtools index Control.bam
```

You can download hg38 bowtie2 index files:

```
mkdir hg38
cd hg38
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/hg38_bowtie2.tar.gz
tar -xzvf hg38_bowtie2.tar.gz
cd ..
```
## 5. (Optional) Correct GC content bias

deepTools can be used to correct GC bias on bam file.
```
(i)	computeGCBias -b [IN_BAM] –g [GENOME] --effectiveGenomeSize [GENOME_SIZE] –GCbiasFrequenciesFile [TEXT_OUT] [-p [THREADS]]
(ii)	correctGCBias -b [IN_BAM] –o [OUT_BAM] –g [GENOME] --effectiveGenomeSize [GENOME_SIZE]–GCbiasFrequenciesFile [TEXT_OUT] [-p [THREADS]]
```
## 6. Count UV radiation induced DNA damage

This step completes raw data processing and generates important output files that will be required in multiple steps in the following workflow. 

```
usage: cpdseqer bam2dinucleotide [-h] -i [INPUT] -g [FASTA] [-q [MAPPING_QUALITY]] [-m [MIN_COVERAGE]] [-u] [-t] -o [OUTPUT]

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
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/UV.bam
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/UV.bam.bai
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/Control.bam
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/Control.bam.bai
```
## 7. (Optional) Generate text files to inform estimated lesions along chromosomes
(i)	Generate a binwise dinucleotide site summary based on the reference genome
```
cpdseqer fasta2bincount -i [INPUT_FA] [-b [BLOCK]] -o [OUTPUT]
```
(ii)	Generate a binwise dinucleotide read-count summary based on CPD read count information
```
cpdseqer dinucleotide2bincount -i [INPUT_DI] [-g [GENOME]] [-b [BLOCK]] -o [OUTPUT] 
```
## 8. (Optional) Subtract out short tandem repeat regions or narrow down to the interested gnomic regions 
``` 
cpdseqer filter –i [INPUT] –c [COORDINATE] –o [OUTPUT_PREFIX] [-m {subtract,intersect}]
```
## 9. (Optional) Estimate sample-wise normalization factors
```
cpdseqer size_factor -i [INPUT] -o [OUTPUT_PREFIX] [--calc_type {site_union,chrom_dinucleotide} }]
```
## 10. Quality control
QC based on dinucleotide count results can be performed using the following command:
```
cpdseqer qc [-h] -i [INPUT] -o [OUTPUT] [-n [NAME]] [--count_type {rCnt,sCnt}}] [-g [GENOME]] [-s [SIZE_FACTOR_FILE]]
```
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

The [dinucleotide.list](https://cqsweb.app.vumc.org/Data/cpdseqer/data/dinucleotide.list) contains two columns indicate dinucleotide file and sample name (separated by tab).

|||
|---|---|
|Control.bed.bgz|Control|
|UV.dinucleotide.bed.bgz|UV|

You can download example files as following scripts.
```
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/UV.bed.bgz
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/UV.bed.bgz.tbi
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/UV.count
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/Control.bed.bgz
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/Control.bed.bgz.tbi
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/Control.count
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/dinucleotide.list
```

The [coordinates.list](https://cqsweb.app.vumc.org/Data/cpdseqer/data/coordinates.list) contains two columns indicate coordinate bed file and category name (separated by tab)

|||
|---|---|
|hg38_promoter.bed|Promoter|
|hg38_tf.bed|TFBinding|

You can download example files as following scripts. The hg38_promoter.bed contains three columns only. So, Promoter (from  coordinates.list definition) will be used as category name for all entries in the hg38_promoter.bed. The hg38_tf.bed contians four columns. The forth column in hg38_tf.bed indicates TF name which will be used as category name (--category_index 3) instead of TFBinding.

```
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/coordinates.list
wget https://github.com/shengqh/cpdseqer/raw/master/data/hg38_promoter.bed
wget https://github.com/shengqh/cpdseqer/raw/master/data/hg38_tf.bed

```
## 11. Generate a genome-wide UV damage distribution map

```
usage: cpdseqer fig_genome [-h] -i [INPUT] [-b [BLOCK]] [-d [DB]] [-n [{None,Total,LocalGC}]] -o [OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT], --input [INPUT]
                        Input dinucleotide list file, first column is file location, second column is file name
  -b [BLOCK], --block [BLOCK]
                        Block size for summerize dinucleotide count (default 100000)
  -d [DB], --db [DB]    Input database version, hg38 or hg19 (default hg38)
  -n [{None,Total,LocalGC}], --norm_type [{None,Total,LocalGC}]
                        Normalization type
  -o [OUTPUT], --output [OUTPUT]
                        Output file prefix
```

for example:

```
cpdseqer fig_genome -i dinucleotide.list -d hg38 -n Total -o output_prefix
```

The [dinucleotide.list](https://cqsweb.app.vumc.org/Data/cpdseqer/data/dinucleotide.list) contains two columns indicate dinucleotide file and sample name (separated by tab).

## 12. Draw dinucleotide pileup figure in a specific genomic region type

```
cpdseqer uv_comp_genome [-h] -i [INPUT] -o [OUTPUT] [–g [GENOME]] [--count_type [COUNT_TYPE]] [-s [SIZE_FACTOR_FILE]]
```
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

You can download hg19/hg38 and yeast background file from:

```
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/hg19_Naked.bed.bgz
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/hg19_Naked.bed.bgz.tbi
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/hg19_Naked.count
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/hg38_Naked.bed.bgz
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/hg38_Naked.bed.bgz.tbi
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/hg38_Naked.count
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/Yeast_Naked.bed.bgz
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/Yeast_Naked.bed.bgz.tbi
wget https://cqsweb.app.vumc.org/Data/cpdseqer/data/Yeast_Naked.count
```

## 13. Compare UV radiation damage of sample(s) against the reference genome background

```
cpdseqer uv_comp_genome [-h] -i [INPUT] -o [OUTPUT] [–g [GENOME]] [--count_type [COUNT_TYPE]] [-s [SIZE_FACTOR_FILE]]
```

## 14. Compare UV damage of sample(s) against reference genome background within a  specific region type

```
cpdseqer uv_comp_genome_region [-h] -i [INPUT] -o [OUTPUT] -c [COORDINATE_FILE]  -f [FASTA] [--add_chr] [--space] [--count_type [COUNT_TYPE]] [-s [SIZE_FACTOR_FILE]]
```

## 15. Compare UV damage between two regions for one or multiple samples

```
cpdseqer uv_comp_regions [-h] -i [INPUT] -o [OUTPUT] -c1 [COORDINATE_FILE1] -c2 [COORDINATE_FILE2] -f [FASTA] [--add_chr] [--space] [--count_type [COUNT_TYPE]]
```

## 16. Compare genome-wide UV damage between two groups of samples

```
cpdseqer uv_comp_groups [-h] -i1 [INPUT1] -i2 [INPUT2] -o [OUTPUT] [--count_type [COUNT_TYPE]] [-s [SIZE_FACTOR_FILE]]
```

## 17. Compare UV damage between two groups of samples within a specific region type

```
cpdseqer uv_comp_groups_region [-h] -i1 [INPUT1] -i2 [INPUT2] -o [OUTPUT] -c [COORDINATE_FILE] [--add_chr] [--space] [--count_type [COUNT_TYPE]] [-s [SIZE_FACTOR_FILE]]
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
wget -r --no-parent --reject "index.html*" https://cqsweb.app.vumc.org/Data/cpdseqer/
cd cqsweb.app.vumc.org/download1/cpdseqer/data
tar -xzvf GRCh38.p13.bowtie2.tar.gz
rm GRCh38.p13.bowtie2.tar.gz
cd ..
unzip test.zip
rm test.zip
```
