import gzip
import os
import os.path
import sys
import shutil
from collections import OrderedDict

from .common_utils import check_file_exists, get_reference_start, runCmd, readFileMap, checkFileMap
 
def demultiplex(logger, inputFile, outputFolder, configFile, args):
  check_file_exists(inputFile)

  logger.info("reading barcode file: " + configFile + " ...")
  sampleSeqMap = readFileMap(configFile)
  seqSampleMap = {sampleSeqMap[k]:k for k in sampleSeqMap.keys()}
  
  seqFileMap = {}
  barcodeLength = 0
  for seq in seqSampleMap.keys():
    barcodeLength = len(seq)
    seqFile = outputFolder + "/" + seqSampleMap[seq] + ".fastq.gz"
    seqFileMap[seq] = gzip.open(seqFile, 'wt')
  
  logger.info("reading input file: " + inputFile + " ...")
  count = 0
  
  if inputFile.endswith(".gz"):
    fin = gzip.open(inputFile, 'rt')
  else:
    fin = open(inputFile, "rt")

  with fin:
    while(True):
      query = fin.readline()
      if not query:
        break
      
      seq = fin.readline().rstrip()
      skipline = fin.readline()
      score = fin.readline().rstrip()
      
      count = count + 1
      if count % 100000 == 0:
        logger.info("%d processed" % count)

      barcode = seq[0:barcodeLength]
      if barcode in seqFileMap:
        newseq = seq[barcodeLength:-1]
        newscore = score[barcodeLength:-1]
        fout = seqFileMap[barcode]
        fout.write(query)
        fout.write(newseq + '\n')
        fout.write(skipline)
        fout.write(newscore + '\n')
        
  for seq in seqFileMap.keys():
    seqFileMap[seq].close()
  
  logger.info("demultiplex done.")
