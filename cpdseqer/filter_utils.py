import argparse
import sys
import logging
import os
import gzip
import math
import subprocess
from Bio import bgzf
from collections import OrderedDict
from _ctypes import ArgumentError

from .common_utils import check_file_exists, runCmd, write_count_file, check_tool_exists, dinucleotide_to_count

def filter(logger, source_file, target_file, output_prefix, method):
  check_file_exists(source_file)
  check_file_exists(target_file)

  check_tool_exists("bedtools")
  check_tool_exists("tabix")

  tmp_bed = output_prefix + ".tmp.bed"

  method_command = "bedtools %s -a \"%s\" -b \"%s\" > %s" % (method, source_file, target_file, tmp_bed)
  runCmd(method_command, logger)

  if not os.path.isfile(tmp_bed):
    raise Exception("bedtools failed, no output file generated.")

  tmp_file = output_prefix + ".tmp.bed.bgz"
  logger.info("Writing dinucleotide to " + tmp_file + " ...")
  with bgzf.BgzfWriter(tmp_file, "wb") as fout:
    with open(tmp_bed, "rt") as fin:
      for line in fin:
        fout.write(line)
  os.remove(tmp_bed)

  output_file = output_prefix + ".bed.bgz"
  if os.path.exists(output_file):
    os.remove(output_file)
  os.rename(tmp_file, output_file)
  runCmd("tabix -p bed %s " % output_file, logger)

  count_file = output_prefix + ".count"
  dinucleotide_to_count(logger, output_file, count_file)

  logger.info("done.")
