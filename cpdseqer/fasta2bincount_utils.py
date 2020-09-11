import pysam
import os
import os.path
import errno
import sys
import shutil
from collections import OrderedDict
from Bio import SeqIO, bgzf
from Bio.Seq import Seq

from .common_utils import MUT_LEVELS, DINU_LEVELS, check_file_exists, runCmd, read_coordinate_file, dinucleotide_to_count

def fasta2bincount(logger, fasta_file, output_file, block):
  check_file_exists(fasta_file)

  dinus = [di for di in MUT_LEVELS]
  dinus.extend([di for di in DINU_LEVELS if di not in MUT_LEVELS])

  tmp_file = output_file + ".tmp"
  with open(tmp_file, "wt") as fout:
    fout.write("#chrom\tstart\tend\t%s\n" % "\t".join(dinus))
    with open(fasta_file, "rt") as fin:  
      for record in SeqIO.parse(fin,'fasta'):
        id = record.id

        logger.info("Counting dinucleotide of chromosome " + id + " ...")

        seq = str(record.seq)
        istart = 0
        while(istart < len(seq)):
          iend = min(istart + block, len(seq)) - 2
          count_dic = {di:0 for di in dinus}
          for si in range(istart, iend):
            dinu = seq[si:(si+2)].upper()
            if dinu in count_dic:
              count_dic[dinu] += 1

          fout.write("%s\t%d\t%d\t%s\n" % (id, istart, iend+2, "\t".join([str(count_dic[dinu]) for dinu in dinus])))
          istart += block

  if os.path.exists(output_file):
    os.remove(output_file)
  os.rename(tmp_file, output_file)

  logger.info("done.")
