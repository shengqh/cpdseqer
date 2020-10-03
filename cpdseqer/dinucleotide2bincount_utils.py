import pysam
import os
import os.path
import errno
import sys
import shutil
import tabix
from collections import OrderedDict
from Bio import SeqIO, bgzf
from Bio.Seq import Seq
from .CategoryItem import CategoryItem

from .common_utils import MUT_LEVELS, DINU_LEVELS, check_file_exists, read_chromosomes, read_chromosome_length, is_major_chromosome_without_chrM, get_blocks

def dinucleotide2bincount(logger, dinucleotide_file, output_file, block, genome):
  idxFile = dinucleotide_file + ".tbi"
  countFile = dinucleotide_file.replace(".bed.bgz", ".count")

  check_file_exists(dinucleotide_file)
  check_file_exists(idxFile)
  check_file_exists(countFile)

  if os.path.isfile(genome):
    chromInfo_file = genome
  elif genome == "hg19":
    chromInfo_file = os.path.join(os.path.dirname(__file__), "data", "chromInfo_hg19.txt")
  elif genome == 'hg38':
    chromInfo_file = os.path.join(os.path.dirname(__file__), "data", "chromInfo_hg38.txt")
  else:
    raise Exception("I don't understand genome, it can be hg19/hg38 or chromsome length file: %s" % genome)

  logger.info("Reading chromosome length from: %s ..." % chromInfo_file)
  
  chrom_length_map = read_chromosome_length(chromInfo_file, is_major_chromosome_without_chrM)
  cat_items = get_blocks(chrom_length_map, block)

  logger.info("Processing %s ..." % dinucleotide_file)

  dinus = [di for di in MUT_LEVELS]
  dinus.extend([di for di in DINU_LEVELS if di not in MUT_LEVELS])

  tmp_file = output_file + ".tmp"
  with open(tmp_file, "wt") as fout:
    fout.write("#chrom\tstart\tend\t%s\t%s\n" % ("\t".join("ReadCount_" + dinu for dinu in dinus), "\t".join("SiteCount_" + dinu for dinu in dinus)))

    tb = tabix.open(dinucleotide_file)

    step = 100000000 / block
    if step < 0:
      step = 1

    count = 0
    for cat_item in cat_items:
      count = count + 1
      if count % step == 0:
        logger.info("%d / %d" % (count, len(cat_items)))

      #logger.info("Processing %s:%d-%d ..." %(cat_item.reference_name, cat_item.reference_start, cat_item.reference_end))
      tbiter = tb.query(cat_item.reference_name, cat_item.reference_start, cat_item.reference_end)
      records = [record for record in tbiter]
      
      total_read_count = {lm:0 for lm in dinus}
      total_site_count = {lm:0 for lm in dinus}
      
      for record in records:
        dinucleotide = record[3]

        diCount = int(record[4])
        if diCount == 0:
          diCount = 1

        total_read_count[dinucleotide] += diCount
        total_site_count[dinucleotide] += 1

      fout.write("%s\t%d\t%d\t%s\t%s\n" % (cat_item.reference_name, cat_item.reference_start, cat_item.reference_end, "\t".join(str(total_read_count[dinu]) for dinu in dinus), "\t".join(str(total_site_count[dinu]) for dinu in dinus)))

  if os.path.exists(output_file):
    os.remove(output_file)
  os.rename(tmp_file, output_file)

  logger.info("done.")
