import argparse
import sys
import logging
import os
import gzip
import math
import subprocess
from _ctypes import ArgumentError

def bowtie2(logger, inputFile, outputFile, databasePrefix, thread):
  logger.info("Start bowtie2 ...")

  if not os.path.exists(databasePrefix + ".rev.2.bt2") :
    raise ArgumentError("Bowtie2 index not exists: %s" % databasePrefix)

  logfile = outputFile + ".log"

  tmpFile = outputFile + ".sam.tmp"
  if not os.path.isfile(tmpFile):
    with open(tmpFile, "wt") as fout:
      with open(logfile, "wt") as flog:
        subprocess.call(['bowtie2', '-p', str(thread), '-x', databasePrefix, '-U', inputFile, '-S', tmpFile], 
          stderr=flog)

  if not os.path.isfile(tmpFile):
    raise Exception("Bowtie2 failed, no output file generated, check log file: %s" % logfile)

  samFile = outputFile + ".sam"
  if not uniqueOnly:
    os.rename(tmpFile, samFile)

  #tobe finish if we need it.
  logger.info("done")
