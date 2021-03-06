import pysam
import os
import os.path
import errno
import sys
import shutil
from collections import OrderedDict
from Bio import SeqIO, bgzf
from Bio.Seq import Seq

from .DinucleotideItem import DinucleotideItem

from .common_utils import check_file_exists, get_reference_start, runCmd, readFileMap, checkFileMap, dinucleotide_to_count

def bam2dinucleotide(logger, bamFile, output_prefix, genomeFastaFile, mappingQuality=20, uniqueOnly=False, minCoverage=1, isTest=False):
  check_file_exists(bamFile)
  check_file_exists(genomeFastaFile)

  logger.info("reading bam file %s ..." % bamFile )
  dinuItems = []
  count = 0
  with pysam.AlignmentFile(bamFile, "rb") as sf:
    for s in sf.fetch():
      count = count + 1
      if count % 1000000 == 0:
        logger.info(count)

        if isTest:
          break
        
      if s.is_unmapped:
        continue

      if s.is_paired and (not s.is_read1):
        continue

      if s.mapping_quality < mappingQuality:
        continue
        
      if uniqueOnly:
        isUnique=True
        for tag in s.tags:
          if tag[0] == 'XS':
            isUnique=False
            break
        if not isUnique:
          continue

      if s.is_reverse:
        dinuItems.append(DinucleotideItem(s.reference_name, s.reference_end, s.reference_end + 2, s.query_name, s.mapping_quality, "-", "" ))
      else:
        dinuItems.append(DinucleotideItem(s.reference_name, s.reference_start - 2, s.reference_start, s.query_name, s.mapping_quality, "+", "" ))
  
  chrDinuMap = OrderedDict()
  for di in dinuItems:
    chrDinuMap.setdefault(di.reference_name, []).append(di)
  
  for chr in chrDinuMap.keys():
    values = chrDinuMap[chr]
    logger.info("sort %d dinucleotides in chromosome %s..." % (len(values), chr ))
    values.sort(key=get_reference_start)
    logger.info("combine %d dinucleotides in chromosome %s..." % (len(values), chr ) )
    idx = len(values) - 1
    deleteList = set()
    while(idx > 0):
      curDinu = values[idx]
      prev = idx -1
      while(prev >= 0):
        prevDinu = values[prev]
        if curDinu.reference_start != prevDinu.reference_start:
          break
        if curDinu.strand == prevDinu.strand:
          prevDinu.count = prevDinu.count + curDinu.count
          deleteList.add(idx)
          break
        prev = prev - 1
      idx = idx -1
    chrDinuMap[chr] = [i for j, i in enumerate(values) if j not in deleteList]
    logger.info("after combine, there is %d dinucleotides in chromosome %s..." % (len(chrDinuMap[chr]), chr ) )

  if minCoverage > 1:
    for chr in chrDinuMap.keys():
      values = chrDinuMap[chr]
      chrDinuMap[chr] = [v for v in values if v.count >= minCoverage]

  with open(genomeFastaFile, "rt") as fin:  
    for record in SeqIO.parse(fin,'fasta'):
      id = record.id
      logger.info("Filling dinucleotide of " + id + " ...")

      if id in chrDinuMap.keys():
        seq = str(record.seq)
        seqlen = len(seq)
        chrDinuItems = chrDinuMap[id]
        for di in chrDinuItems:
          if di.reference_start >= 0 and di.reference_end <= seqlen:
            dinu = seq[di.reference_start:di.reference_end].upper()
            if di.strand == "+":
              dinu = str(Seq(dinu).reverse_complement())
            di.dinucleotide = dinu
  
  tmp_file = output_prefix + ".tmp.bed.bgz"
  logger.info("Writing dinucleotide to " + tmp_file + " ...")
  with bgzf.BgzfWriter(tmp_file, "wb") as fout:
    for chrom in chrDinuMap.keys():
      diList = chrDinuMap[chrom]
      for s in diList:
        if (s.dinucleotide != "") and (not 'N' in s.dinucleotide):
          fout.write("%s\t%d\t%d\t%s\t%d\t%s\n" % (s.reference_name, s.reference_start, s.reference_end, s.dinucleotide, s.count, s.strand))
  
  output_file = output_prefix + ".bed.bgz"
  if os.path.exists(output_file):
    os.remove(output_file)
  os.rename(tmp_file, output_file)
  runCmd("tabix -p bed %s " % output_file, logger)

  count_file = output_prefix + ".count"
  dinucleotide_to_count(logger, output_file, count_file)

  logger.info("done.")
