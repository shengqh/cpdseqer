# CPDSeqer

This package is used to do CPD sequence data analysis.

# Table of Contents
1. [Prerequisites](#Prerequisites)
2. [Installation](#Installation)
3. [Step by step tutorial](#step-by-step-tutorial)
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
wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
tar -jxvf samtools-1.10.tar.bz2
cd samtools-1.10
./configure
make
sudo make install
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

# Step by step tutorial

We provide example data and code for test. It requires 34G for downloading and additional 20G for testing. Remember to come back to folder cpdseqer after each step.

```
wget -r --no-parent --reject "index.html*" --reject-regex "bed_files" https://cqsweb.app.vumc.org/Data/cpdseqer/
cd cqsweb.app.vumc.org/Data/cpdseqer/
```

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

For example:

```
cd T01_demultiplex/
cpdseqer demultiplex -i ../data/example.fastq.gz -o . -b barcode.txt
```

The [barcode.txt](https://cqsweb.app.vumc.org/Data/cpdseqer/data/barcode.txt) contains two columns indicate barcode and sample name (separated by tab).

|||
|---|---|
|ATCGCGAT|Control|
|GAACTGAT|UV|

An alternative method to perform demultiplex is using [Je](https://gbcs.embl.de/portal/tiki-index.php?page=Je).

```
cd T01_demultiplex/je
je demultiplex F1=../../data/example.fastq.gz BF=barcode_je.txt O=.
```

The [JE_BARCODEFILE](https://cqsweb.app.vumc.org/Data/cpdseqer/T01_demultiplex/je/barcode_je.txt) contains three columns indicate sample name, barcode and sample filename (separated by tab).

||||
|---|---|---|
|Control|ATCGCGAT|Control.fastq.gz|
|UV|GAACTGAT|UV.fastq.gz|

## 2. General QC
Without going into detail, we recommend previously established methods such as FASTQC and QC3 for this general QC step.

```
cd T02_fastqc/
fastqc -o . -t 1 --extract ../data/Control.fastq.gz
```

## 3. Generate index files for reference genome (FASTA) files

```
cd T03_build_bowtie2_index
bowtie2-build --threads 8 ../data/bowtie2_index_2.3.5.1/GRCh38.primary_assembly.genome.fa GRCh38.primary_assembly.genome
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

For example,

```
cd T04_bowtie2
bowtie2 -p 8 -x ../data/bowtie2_index_2.3.5.1/GRCh38.primary_assembly.genome -U ../data/Control.fastq.gz | samtools sort -@ 8 -m 4G -o Control.bam -
samtools index Control.bam
```

## 5. (Optional) Correct GC content bias

deepTools can be used to correct GC bias on bam file.
```
(i)	computeGCBias -b [IN_BAM] –g [GENOME] --effectiveGenomeSize [GENOME_SIZE] –GCbiasFrequenciesFile [TEXT_OUT] [-p [THREADS]]
(ii)	correctGCBias -b [IN_BAM] –o [OUT_BAM] –g [GENOME] --effectiveGenomeSize [GENOME_SIZE]–GCbiasFrequenciesFile [TEXT_OUT] [-p [THREADS]]
```

For example,

```
cd T05_Correct_GC/T05_01_calc
computeGCBias -p 8 -b ../../data/Control.bam --effectiveGenomeSize 2913022398 -g ../../data/hg38.2bit -l 200 --GCbiasFrequenciesFile control_freq.txt --biasPlot control.png

cd ../T05_02_correct/
correctGCBias -p 8 -b ../../data/Control.bam --effectiveGenomeSize 2913022398 -g ../../data/hg38.2bit --GCbiasFrequenciesFile ../T05_01_calc/control_freq.txt -o control.corrected.bam
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

For example:

```
cd T06_bam2dinucleotide

cpdseqer bam2dinucleotide -i ../data/Control.bam -g ../data/bowtie2_index_2.3.5.1/GRCh38.primary_assembly.genome.fa -o Control
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

For example,

```
cd T07_bincount/T07_01_fasta2bincount/
cpdseqer fasta2bincount -i ../../data/bowtie2_index_2.3.5.1/GRCh38.primary_assembly.genome.fa -b 100000 -o GRCh38.b100000.txt
```

Once you come back to cpdseqer folder,

```
cd T07_bincount/T07_02_dinucleotide2bincount/
cpdseqer dinucleotide2bincount -i ../../data/Control.bed.bgz -g hg38 -b 100000 -o Control.b100000.txt
```

## 8. (Optional) Subtract out short tandem repeat regions or narrow down to the interested gnomic regions 

``` 
cpdseqer filter –i [INPUT] –c [COORDINATE] –o [OUTPUT_PREFIX] [-m {subtract,intersect}]
```

For example,

```
cd T08_filter
cpdseqer filter -i ../data/Control.bed.bgz -c ../data/hg38_UTR3.bed.gz -m intersect -o Control.UTR3
cpdseqer filter -i ../data/Control.bed.bgz -c ../data/hg38_UTR3.bed.gz -m subtract -o Control.noUTR3
```

## 9. (Optional) Estimate sample-wise normalization factors
```
cpdseqer size_factor -i [INPUT] -o [OUTPUT_PREFIX] [--calc_type {site_union,chrom_dinucleotide} }]
```

For example,

```
cd T09_sizefactor
cpdseqer size_factor -i dinucleotide.list --calc_type chrom_dinucleotide -o sf
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

For example:

```
cd T10_qc
cpdseqer qc -i single_file.list -n single -o output_single
cpdseqer qc -i multi_file.list -n multi -o output_multi
```

## 11. Generate a genome-wide UV damage distribution map

```
usage: cpdseqer fig_genome [-h] -i [INPUT] [-b [BLOCK]]
                           [-n [{None,Total,LocalGC}]] -o [OUTPUT]
                           [-g [GENOME]]

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT], --input [INPUT]
                        Input dinucleotide list file, first column is file
                        location, second column is file name
  -b [BLOCK], --block [BLOCK]
                        Block size for summerize dinucleotide count (default
                        100000)
  -n [{None,Total,LocalGC}], --norm_type [{None,Total,LocalGC}]
                        Normalization type
  -o [OUTPUT], --output [OUTPUT]
                        Output file prefix
  -g [GENOME], --genome [GENOME]
                        Input reference genome, hg38/hg19 (default hg38) or
                        chromosome length file
```

For example:

```
cd T11_fig_genome
cpdseqer fig_genome -i dinucleotide.list -g hg38 -o cpd_genome -n LocalGC -b 1000000
```

The [dinucleotide.list](https://cqsweb.app.vumc.org/Data/cpdseqer/T11_fig_genome/dinucleotide.list) contains two columns indicate dinucleotide file and sample name (separated by tab).

|||
|---|---|
|../data/Control.bed.bgz|Control|
|../data/UV.dinucleotide.bed.bgz|UV|

## 12. Draw dinucleotide pileup figure in a specific genomic region type

```
usage: cpdseqer fig_position [-h] -i [INPUT] -c [COORDINATE_FILE]
                             [-b [BACKGROUND_FILE]] [--space] [--add_chr] [-t]
                             -o [OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT], --input [INPUT]
                        Input dinucleotide list file, first column is file
                        location, second column is file name
  -c [COORDINATE_FILE], --coordinate_file [COORDINATE_FILE]
                        Input coordinate bed file (can use short name
                        hg38/hg19 as default nucleosome file)
  -b [BACKGROUND_FILE], --background_file [BACKGROUND_FILE]
                        Background dinucleotide file
  --space               Use space rather than tab in coordinate files
  --add_chr             Add chr in chromosome name in coordinate file
  -t, --test            Test the first 10000 coordinates only
  -o [OUTPUT], --output [OUTPUT]
                        Output file prefix
```

For example, we will calculate the dinucleotide position in 3'UTR:

```
cd T12_fig_position
cpdseqer fig_position -c ../data/hg38_UTR3.bed.gz -b ../data/hg38_Naked.bed.bgz -i dinucleotide.list -o position
```

Here, you can input absolute coordinate file, or hg38/hg19. hg38 and hg19 indicates the nucleosome coordinate files which can also be downloaded by:

```
wget https://github.com/shengqh/cpdseqer/raw/master/cpdseqer/data/nucleosome_hg19_interval.zip
wget https://github.com/shengqh/cpdseqer/raw/master/cpdseqer/data/nucleosome_hg38_interval.zip

```

## 13. Compare UV radiation damage of sample(s) against the reference genome background

```
usage: cpdseqer uv_comp_genome [-h] -i [INPUT] [--count_type [{rCnt,sCnt}]] -o
                               [OUTPUT] [-g [GENOME]] [-s [SIZE_FACTOR_FILE]]

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT], --input [INPUT]
                        Input count list file, first column is file location,
                        second column is file name
  --count_type [{rCnt,sCnt}]
                        Input count type, rCnt/sCnt (read count/site count,
                        default rCnt)
  -o [OUTPUT], --output [OUTPUT]
                        Output file prefix
  -g [GENOME], --genome [GENOME]
                        Input reference genome, hg38/hg19/saccer3 (default
                        hg38)
  -s [SIZE_FACTOR_FILE], --size_factor_file [SIZE_FACTOR_FILE]
                        Input size factor file for normalization
```

For example, 

```
cd T13_uv_comp_genome
cpdseqer uv_comp_genome -i count.list -o uv_comp_genome -g hg38 --count_type sCnt -s ../data/chrom_dinucleotide.sizefactor.txt
```

The [count.list](https://cqsweb.app.vumc.org/Data/cpdseqer/T13_uv_comp_genome/count.list) contains two columns indicate dinucleotide count file and sample name (separated by tab).

|||
|---|---|
|../data/Control.count|Control|
|../data/UV.count|UV|

## 14. Compare UV damage of sample(s) against reference genome background within a  specific region type

```
usage: cpdseqer uv_comp_genome_region [-h] -i [INPUT] [-c [COORDINATE_FILE]]
                                      [--space] [--add_chr] [-f [FASTA]]
                                      [--count_type [{rCnt,sCnt}]] -o [OUTPUT]
                                      [-s [SIZE_FACTOR_FILE]]

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT], --input [INPUT]
                        Input dinucleotide list file, first column is file
                        location, second column is file name
  -c [COORDINATE_FILE], --coordinate_file [COORDINATE_FILE]
                        Input coordinate bed file (can use short name
                        hg38/hg19 as default nucleosome file)
  --space               Use space rather than tab in coordinate file
  --add_chr             Add chr to chromosome name in coordinate file
  -f [FASTA], --fasta [FASTA]
                        Input reference genome fasta file
  --count_type [{rCnt,sCnt}]
                        Input count type, rCnt/sCnt (read count/site count,
                        default rCnt)
  -o [OUTPUT], --output [OUTPUT]
                        Output file prefix
  -s [SIZE_FACTOR_FILE], --size_factor_file [SIZE_FACTOR_FILE]
                        Input size factor file for normalization
```

For example,

```
cd T14_uv_comp_genome_region
cpdseqer uv_comp_genome_region \
  -i dinucleotide.list \
  -f ../data/bowtie2_index_2.3.5.1/GRCh38.primary_assembly.genome.fa \
  -c ../data/nucleosome_hg38_interval.txt \
  -o uv_comp_genome_region \
  --space  \
  -s ../data/chrom_dinucleotide.sizefactor.txt
```

## 15. Compare UV damage between two regions for one or multiple samples

```
usage: cpdseqer uv_comp_regions [-h] -i [INPUT] [-c1 [COORDINATE_FILE1]]
                                [-c2 [COORDINATE_FILE2]] [--space] [--add_chr]
                                [-f [FASTA]] [--count_type [{rCnt,sCnt}]] -o
                                [OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT], --input [INPUT]
                        Input dinucleotide list file, first column is file
                        location, second column is file name
  -c1 [COORDINATE_FILE1], --coordinate_file1 [COORDINATE_FILE1]
                        Input coordinate bed file 1
  -c2 [COORDINATE_FILE2], --coordinate_file2 [COORDINATE_FILE2]
                        Input coordinate bed file 2
  --space               Use space rather than tab in coordinate file
  --add_chr             Add chr to chromosome name in coordinate file
  -f [FASTA], --fasta [FASTA]
                        Input reference genome fasta file
  --count_type [{rCnt,sCnt}]
                        Input count type, rCnt/sCnt (read count/site count,
                        default rCnt)
  -o [OUTPUT], --output [OUTPUT]
                        Output file prefix
```

For example,

```
cd T15_uv_comp_regions
cpdseqer uv_comp_regions \
  -i dinucleotide.list \
  -f ../data/bowtie2_index_2.3.5.1/GRCh38.primary_assembly.genome.fa \
  -c1 hg38_promoter.bed \
  -c2 hg38_tf.bed \
  -o uv_comp_regions
```

## 16. Compare genome-wide UV damage between two groups of samples

```
usage: cpdseqer uv_comp_groups [-h] -i1 [INPUT1] -i2 [INPUT2]
                               [--count_type [{rCnt,sCnt}]] -o [OUTPUT]
                               [-s [SIZE_FACTOR_FILE]]

optional arguments:
  -h, --help            show this help message and exit
  -i1 [INPUT1], --input1 [INPUT1]
                        Input CPD count list file 1, first column is file
                        location, second column is file name
  -i2 [INPUT2], --input2 [INPUT2]
                        Input CPD count list file 2, first column is file
                        location, second column is file name
  --count_type [{rCnt,sCnt}]
                        Input count type, rCnt/sCnt (read count/site count,
                        default rCnt)
  -o [OUTPUT], --output [OUTPUT]
                        Output file prefix
  -s [SIZE_FACTOR_FILE], --size_factor_file [SIZE_FACTOR_FILE]
                        Input size factor file for normalization
```

For example,

```
cd T16_uv_comp_groups
cpdseqer uv_comp_groups \
  -i1 control.list \
  -i2 case.list \
  -o uv_comp_groups \
  --count_type sCnt \
  -s ../data/chrom_dinucleotide.sizefactor.txt
```

## 17. Compare UV damage between two groups of samples within a specific region type

```
cpdseqer uv_comp_groups_region [-h] -i1 [INPUT1] -i2 [INPUT2] -o [OUTPUT] -c [COORDINATE_FILE] [--add_chr] [--space] [--count_type [COUNT_TYPE]] [-s [SIZE_FACTOR_FILE]]
```

For example,

```
cd T17_uv_comp_groups_region
cpdseqer uv_comp_groups_region \
  -i1 control.list \
  -i2 case.list \
  -c ../data/nucleosome_hg38_interval.txt \
  -o uv_comp_groups_region \
  --count_type sCnt \
  --space \
  -s ../data/chrom_dinucleotide.sizefactor.txt
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
