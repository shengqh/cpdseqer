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

from .common_utils import check_file_exists, runCmd, write_count_file, check_tool_exists

def filter(logger, sourceFile, targetFile, outputPrefix, method):
  check_file_exists(sourceFile)
  check_file_exists(targetFile)

  check_tool_exists("bedtools")
  check_tool_exists("tabix")

  tmpFile = outputPrefix + ".tmp.bed"

  method_command = "bedtools %s -a \"%s\" -b \"%s\" > %s" % (method, sourceFile, targetFile, tmpFile)
  runCmd(method_command, logger)

  if not os.path.isfile(tmpFile):
    raise Exception("bedtools failed, no output file generated.")

  outputFile = outputPrefix + ".bed.bgz"
  countMap = OrderedDict()
  logger.info("Writing dinucleotide to " + outputFile + " ...")
  with bgzf.BgzfWriter(outputFile, "wb") as fout:
    with open(tmpFile, "rt") as fin:
      for line in fin:
        fout.write(line)
        parts = line.split('\t')
        chrom = parts[0]
        dinucleotide = parts[3]
        count = int(parts[4])
        chromMap = countMap.setdefault(chrom, {})
        countVec = chromMap.setdefault(dinucleotide, [0,0])
        countVec[0] = countVec[0] + count
        countVec[1] = countVec[1] + 1
  
  countFile = outputPrefix + ".count"
  write_count_file(logger, countFile, countMap)

  runCmd("tabix -p bed %s " % outputFile, logger)

  os.remove(tmpFile)

  logger.info("done.")
