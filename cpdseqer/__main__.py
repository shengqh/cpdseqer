import argparse
import sys
import os
from pathlib import Path

from .__version__ import __version__
from .common_utils import runCmd, initialize_logger
from .background_utils import background
from .demultiplex_utils import demultiplex
# from .bowtie2_utils import bowtie2
from .bam2dinucleotide_utils import bam2dinucleotide
from .statistic_utils import statistic
from .report_utils import report
from .fig_genome_utils import fig_genome
from .fig_position_utils import fig_position
from .uv_comp_utils import uv_comp_genome, uv_comp_genome_region
from .qc_utils import qc

# CPD-seq data analysis
# CPD-seq sequencing reads were trimmed of barcode sequences and the 3' nucleotide of the sequencing read, 
# and then aligned to the hg19 human genome using the bowtie 2 software44. The resulting alignment files were 
# processed with SAMtools45 and BEDtools46, and custom Perl scripts were used to identify dinucleotide sequence 
# immediately upstream of the 5' end of each sequencing read. The dinucleotide sequence on the opposite strand 
# was extracted as a putative CPD lesion. Background reads associated with non-dipyrimidine sequences, which were 
# likely due to incomplete 3' DNA end blocking or non-specific DNA cleavage by T4 endonuclease V/APE1, were excluded 
# from subsequent analyses. Both positions in the dipyrimidine nucleotide were counted as lesion sites. Three independent 
# CPD-seq experiments mapped CPD lesions in UV-irradiated NHF1 cells (UV cells) and two independent CPD-seq experiments 
# mapped lesions in isolated NHF1 genomic DNA that was UV-irradiated in vitro (UV naked DNA). These biological replicates 
# were combined for most of the analyses. Additionally, in some cases only mutagenic CPD (mCPDs), which are CPD reads 
# associated with TC, CT, or CC dinucleotides, were analyzed.


def main():
  parser = argparse.ArgumentParser(description="CPDseq analysis " + __version__,
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  subparsers = parser.add_subparsers(dest="command")

  # create the parser for the "background" command
  parser_b = subparsers.add_parser('background', help='Calculate background dinucleotide count in genome sequence')
  parser_b.add_argument("-i", "--input", help = "Input fasta file", required=True)
  parser_b.add_argument("-b", "--bed", help = "Input bed file (optional)", default='NONE')
  parser_b.add_argument("-o", "--output", help = "Output file", required=False)

  # create the parser for the "demultiplex" command
  parser_d = subparsers.add_parser('demultiplex', help='Perform demultiplex on raw fastq file')
  parser_d.add_argument('-i', '--input', action='store', nargs='?', help="Input fastq file (gzipped supported)", required=NOT_DEBUG)
  parser_d.add_argument('-o', '--output', action='store', nargs='?', help="Output folder", required=NOT_DEBUG)
  parser_d.add_argument('-b', '--barcode_file', action='store', nargs='?', help='Tab-delimited file, first column is barcode, second column is sample name', required=NOT_DEBUG)

  # create the parser for the "bowtie2" command
  # parser_b = subparsers.add_parser('bowtie2')
  # parser_b.add_argument('-i', '--input', action='store', nargs='?', help='Input single-end fastq/fasta file', required=NOT_DEBUG)
  # parser_b.add_argument('-d', '--database_prefix', action='store', nargs='?', help='Input bowtie2 index prefix', required=NOT_DEBUG)
  # parser_b.add_argument('-t', '--thread', action='store', nargs='?', type=int, default=8, help="Thread number")
  # parser_b.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output file", required=NOT_DEBUG)

  # create the parser for the "bam2dinucleotide" command
  parser_p = subparsers.add_parser('bam2dinucleotide', help='Extract dinucleotide from bam file')
  parser_p.add_argument('-i', '--input', action='store', nargs='?', help='Input BAM file', required=NOT_DEBUG)
  parser_p.add_argument('-g', '--genome_seq_file', action='store', nargs='?', help='Input genome seq file', required=NOT_DEBUG)
  parser_p.add_argument('-q', '--mapping_quality', type=int, default=20, nargs='?', help='Minimum mapping quality of read (default 20)', required=False)
  parser_p.add_argument('-m', '--min_coverage', type=int, default=1, help='The minimum coverage of dinucleotide for counting (default 1)')
  parser_p.add_argument('-u', '--unique_only', action='store_true', help='Use uniquely mapped read only')
  parser_p.add_argument('-t', '--test', action='store_true', help='Test the first 1000000 reads only')
  parser_p.add_argument('-o', '--output', action='store', nargs='?', help="Output file prefix", required=NOT_DEBUG)

  parser_q = subparsers.add_parser('qc', help='Quality control based on dinucleotide count result')
  parser_q.add_argument('-i', '--input', action='store', nargs='?', help='Input dinucleotide prefix', required=NOT_DEBUG)
  parser_q.add_argument('-n', '--name', action='store', nargs='?', help='Input sample name')
  parser_q.add_argument('--count_type', action='store', nargs='?', default="rCnt", help='Input count type, rCnt/sCnt (read count/site count, default rCnt)')
  parser_q.add_argument('-o', '--output', action='store', nargs='?', help="Output file prefix", required=NOT_DEBUG)

  # create the parser for the "statistic" command
  # parser_s = subparsers.add_parser('statistic')
  # parser_s.add_argument('-i', '--input', action='store', nargs='?', help='Input dinucleotide list file, first column is file location, second column is file name', required=NOT_DEBUG)
  # parser_s.add_argument('-c', '--coordinate_list_file', action='store', nargs='?', help='Input coordinate list file, first column is file location, second column is file name', required=NOT_DEBUG)
  # parser_s.add_argument('-s', '--space', action='store_true', help='Use space rather than tab in coordinate files')
  # parser_s.add_argument('--category_index', type=int, default=-1, nargs='?', help='Zero-based category column index in coordinate file')
  # parser_s.add_argument('--add_chr', action='store_true', help='Add chr to chromosome name in coordinate file')
  # parser_s.add_argument('-o', '--output', action='store', nargs='?', help="Output file name", required=NOT_DEBUG)

  # create the parser for the "report" command
  parser_r = subparsers.add_parser('report')
  parser_r.add_argument('-i', '--input', action='store', nargs='?', help='Input dinucleotide list file, first column is file location, second column is file name', required=NOT_DEBUG)
  parser_r.add_argument('-g', '--group', action='store', nargs='?', help='Input group list file, first column is group (0 for control and 1 for case), second column is file name', required=NOT_DEBUG)
  parser_r.add_argument('-b', '--block', type=int, default=100000, nargs='?', help='Block size for summerize dinucleotide count')
  parser_r.add_argument('-d', '--db', action='store', nargs='?', default="hg38", help='Input database version, hg38 or hg19 (default hg38)')
  parser_r.add_argument('-o', '--output', action='store', nargs='?', help="Output file name", required=NOT_DEBUG)

  # create the parser for the "fig_genome" command
  parser_sg = subparsers.add_parser('fig_genome', help="Generate statistic figure on genome/chromosome level")
  parser_sg.add_argument('-i', '--input', action='store', nargs='?', help='Input configuration file', required=NOT_DEBUG)
  parser_sg.add_argument('-b', '--block', type=int, default=100000, nargs='?', help='Block size for summerize dinucleotide count (default 100000)')
  parser_sg.add_argument('-d', '--db', action='store', nargs='?', default="hg38", help='Input database version, hg38 or hg19 (default hg38)')
  parser_sg.add_argument('-o', '--output', action='store', nargs='?', help="Output file name", required=NOT_DEBUG)

  # create the parser for the "fig_position" command
  parser_f = subparsers.add_parser('fig_position', help='Generate statistic figure on relative position in coordinate file')
  parser_f.add_argument('-i', '--input', action='store', nargs='?', help='Input dinucleotide list file, first column is file location, second column is file name', required=NOT_DEBUG)
  parser_f.add_argument('-c', '--coordinate_file', action='store', nargs='?', help='Input coordinate bed file (can use short name hg38/hg19 as default nucleosome file)', required=NOT_DEBUG)
  parser_f.add_argument('-b', '--background_file', action='store', nargs='?', help='Background dinucleotide file')
  parser_f.add_argument('-s', '--space', action='store_true', help='Use space rather than tab in coordinate files')
  parser_f.add_argument('--add_chr', action='store_true', help='Add chr in chromosome name in coordinate file')
  parser_f.add_argument('-t', '--test', action='store_true', help='Test the first 10000 coordinates only')
  parser_f.add_argument('-o', '--output', action='store', nargs='?', help="Output file name", required=NOT_DEBUG)

  parser_u = subparsers.add_parser('uv_comp_genome', help='Compare UV radiation damage between sample(s) and reference genome background')
  parser_u.add_argument('-i', '--input', action='store', nargs='?', help='Input count list file, first column is file location, second column is file name', required=NOT_DEBUG)
  parser_u.add_argument('-d', '--db', action='store', nargs='?', default="hg38", help='Input reference genome version, hg38/hg19/saccer3 (default hg38)')
  parser_u.add_argument('--count_type', action='store', nargs='?', default="rCnt", help='Input count type, rCnt/sCnt (read count/site count, default rCnt)')
  parser_u.add_argument('-o', '--output', action='store', nargs='?', help="Output file prefix", required=NOT_DEBUG)

  parser_u = subparsers.add_parser('uv_comp_genome_region', help='Compare UV radiation damage between sample(s) and reference genome background in specific region')
  parser_u.add_argument('-i', '--input', action='store', nargs='?', help='Input dinucleotide list file, first column is file location, second column is file name', required=NOT_DEBUG)
  parser_u.add_argument('-c', '--coordinate_file', action='store', nargs='?', help='Input coordinate bed file (can use short name hg38/hg19 as default nucleosome file)', required=False)
  parser_u.add_argument('-s', '--space', action='store_true', help='Use space rather than tab in coordinate file')
  parser_u.add_argument('--add_chr', action='store_true', help='Add chr to chromosome name in coordinate file')
  parser_u.add_argument('-d', '--db', action='store', nargs='?', help='Input reference genome fasta file')
  parser_u.add_argument('--count_type', action='store', nargs='?', default="rCnt", help='Input count type, rCnt/sCnt (read count/site count, default rCnt)')
  parser_u.add_argument('-o', '--output', action='store', nargs='?', help="Output file prefix", required=NOT_DEBUG)
  
  if not DEBUG and len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

  args = parser.parse_args()
  print(args)
  
  if args.command == "background":
    logger = initialize_logger(args.output + ".log", args)
    background(args.input, args.bed, args.output)
  elif args.command == "demultiplex":
    if not os.path.isdir(args.output):
      Path(args.output).mkdir(parents=True)
    logger = initialize_logger(args.output + "/cpdseqer_demultiplex.log", args)
    demultiplex(logger, args.input, args.output, args.barcode_file, args)
  # elif args.command == 'bowtie2':
  #   logger = initialize_logger(args.output + ".cpdseqer_bowtie2.log", args)
  #   bowtie2(logger, args.input, args.output, args.database_prefix, args.thread)
  elif args.command == "bam2dinucleotide":
    logger = initialize_logger(args.output + ".log", args)
    bam2dinucleotide(logger, args.input, args.output, args.genome_seq_file, args.mapping_quality, args.unique_only, args.min_coverage, args.test)
  elif args.command == "qc":
    logger = initialize_logger(args.output + ".log", args)
    qc(logger, args.input, args.name, args.output, args.count_type)
  # elif args.command == "statistic":
  #   logger = initialize_logger(args.output + ".log", args)
  #   statistic(logger, args.input, args.output, args.coordinate_list_file, args.space, args.add_chr, args.category_index)
  elif args.command == "report":
    logger = initialize_logger(args.output + ".log", args)
    report(logger, args.input, args.group, args.output, args.block, args.db)
  elif args.command == "fig_genome":
    logger = initialize_logger(args.output + ".log", args)
    fig_genome(logger, args.input, args.output, args.block, args.db)
  elif args.command == "fig_position":
    logger = initialize_logger(args.output + ".log", args)
    fig_position(logger, args.input, args.output, args.coordinate_file, args.background_file, args.space, args.add_chr, args.test)
  elif args.command == "uv_comp_genome":
    logger = initialize_logger(args.output + ".log", args)
    uv_comp_genome(logger, args.input, args.output, args.db, args.count_type)
  elif args.command == "uv_comp_genome_region":
    logger = initialize_logger(args.output + ".log", args)
    uv_comp_genome_region(logger, args.input, args.output, args.db, args.count_type, args.coordinate_file, args.space, args.add_chr)
  
if __name__ == "__main__":
    main()
