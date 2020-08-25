import pysam
import os
import os.path
import errno
import sys
import shutil
from collections import OrderedDict
from Bio import SeqIO, bgzf
from Bio.Seq import Seq

from .common_utils import check_file_exists, runCmd, read_coordinate_file, dinucleotide_to_count

def fasta2dinucleotide(logger, fasta_file, bed_file, output_prefix, is_test=False):
  check_file_exists(fasta_file)
  check_file_exists(bed_file)

  regions = read_coordinate_file(bed_file, "region", checkOverlap=True)
  chromRegionMap = {}
  for region in regions:
    chromRegionMap.setdefault(region.reference_name, []).append(region)

  tmp_file = output_prefix + ".tmp.bed.bgz"
  with bgzf.BgzfWriter(tmp_file, "wb") as fout:
    with open(fasta_file, "rt") as fin:  
      for record in SeqIO.parse(fin,'fasta'):
        id = record.id
        if id not in chromRegionMap:
          continue

        logger.info("Extracting dinucleotide of " + id + " ...")

        seq = str(record.seq)
        catItems = chromRegionMap[id]
        for ci in catItems:
          catSeq = seq[ci.reference_start:ci.reference_end].upper()
          if ci.strand == '-':
            catSeq = str(Seq(catSeq).reverse_complement())
          for si in range(0, len(catSeq) - 2):
            dinu = catSeq[si:(si+2)].upper()
            fout.write("%s\t%d\t%d\t%s\t%d\t%s\n" % (id, ci.reference_start + si, ci.reference_start + si + 2, dinu, 1, ci.strand))

  output_file = output_prefix + ".bed.bgz"
  if os.path.exists(output_file):
    os.remove(output_file)
  os.rename(tmp_file, output_file)
  runCmd("tabix -p bed %s " % output_file, logger)

  count_file = output_prefix + ".count"
  dinucleotide_to_count(logger, output_file, count_file)

  logger.info("done.")
