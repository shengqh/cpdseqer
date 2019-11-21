import argparse
import sys
import logging
import os
from AnalysisUtils import demultiplex, bam2dinucleotide

def initialize_logger(logfile, args):
  logger = logging.getLogger('cpdseqer')
  loglevel = logging.INFO
  logger.setLevel(loglevel)

  formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')    
 
  # create console handler and set level to info
  handler = logging.StreamHandler()
  handler.setLevel(loglevel)
  handler.setFormatter(formatter)
  logger.addHandler(handler)
 
  # create error file handler and set level to error
  handler = logging.FileHandler(logfile, "w")
  handler.setLevel(loglevel)
  handler.setFormatter(formatter)
  logger.addHandler(handler)
 
  return(logger)

def runCommand(command, logger):
  logger.info("run : " + command )
  os.system(command)


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
  parser = argparse.ArgumentParser(description="CPDseq analysis",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  subparsers = parser.add_subparsers(dest="command")

  # create the parser for the "demultiplex" command
  parser_d = subparsers.add_parser('demultiplex')
  parser_d.add_argument('-i', '--input', action='store', nargs='?', help="Input fastq gzipped file", required=NOT_DEBUG)
  parser_d.add_argument('-o', '--output', action='store', nargs='?', help="Output file prefix", required=NOT_DEBUG)
  parser_d.add_argument('-b', '--barcodeFile', action='store', nargs='?', help='Barcode definition file', required=NOT_DEBUG)

  # create the parser for the "bam2dinucleotide" command
  parser_p = subparsers.add_parser('bam2dinucleotide')
  parser_p.add_argument('-i', '--input', action='store', nargs='?', help='Input BAM file', required=NOT_DEBUG)
  parser_p.add_argument('-g', '--genome_seq_file', action='store', nargs='?', help='Input genome seq file', required=NOT_DEBUG)
  parser_p.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output file name", required=NOT_DEBUG)

  # create the parser for the "statistic" command
  parser_s = subparsers.add_parser('statistic')
  parser_s.add_argument('-i', '--input', action='store', nargs='?', help='Input dinucleotide file', required=NOT_DEBUG)
  parser_s.add_argument('-c', '--coordinate_file', action='store', nargs='?', help='Input coordinate file', required=NOT_DEBUG)
  parser_s.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output file name", required=NOT_DEBUG)
  
  if not DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

  args = parser.parse_args()
  
  if args.command == "demultiplex":
    logger = initialize_logger(args.output + "/cpdseqer.log", args)
    demultiplex(logger, args.input, args.output, args.barcodeFile, args)
  elif args.command == "bam2dinucleotide":
    logger = initialize_logger(args.output + ".log", args)
    bam2dinucleotide(logger, args.input, args.output, args.genome_seq_file)
  elif args.command == "statistic":
    logger = initialize_logger(args.output + ".log", args)
    statistic(logger, args.input, args.output, args.coordinate_file)
  
if __name__ == "__main__":
    main()