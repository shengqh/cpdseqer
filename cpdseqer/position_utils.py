import gzip
import pysam
import os
import os.path
import errno
import tabix
import sys
import shutil

from .CategoryItem import CategoryItem
from .common_utils import check_file_exists, get_reference_start, runCmd, readFileMap, checkFileMap, remove_chr

def position(logger, dinucleotideFileList, outputFile, coordinateFile, normalizedByAA=False, useSpace=False, addChr=False):
  check_file_exists(dinucleotideFileList)
  check_file_exists(coordinateFile)

  dinucleotideFileMap = readFileMap(dinucleotideFileList)

  coordinates = []
  delimit = ' ' if useSpace else '\t'
  if (coordinateFile == 'hg38_promoter.bed') or (coordinateFile == 'hg38_tf.bed'):
    if not os.path.exists(coordinateFile):
      coordinateFile = os.path.join(os.path.dirname(__file__), "data", coordinateFile)

  logger.info("Reading coordinate file " + coordinateFile + " ...")
      
  with open(coordinateFile, "rt") as fin:
    for line in fin:
      parts = line.rstrip().split(delimit)
      chrom = "chr" + parts[0] if addChr else parts[0] 
      coordinates.append(CategoryItem(chrom, int(float(parts[1])), int(float(parts[2])), None))

  maxIndex = max([cor.getLength() for cor in coordinates])
        
  finalMap = {dinuName:{index:{} for index in range(0, maxIndex)} for dinuName in dinucleotideFileMap.keys()}
  
  for dinuName in dinucleotideFileMap.keys(): 
    dinucleotideFile = dinucleotideFileMap[dinuName]           
    idxFile = dinucleotideFile + ".tbi"
    countFile = dinucleotideFile + ".count"

    check_file_exists(dinucleotideFile)
    check_file_exists(idxFile)
    check_file_exists(countFile)

    chromMap = {}
    with open(countFile, "rt") as fin:
      fin.readline()
      for line in fin:
        parts = line.split('\t')
        chromMap[parts[0]] = 1

    totalAA = 0
    if normalizedByAA:
      countFile = dinucleotideFile + ".count"
      check_file_exists(countFile)

      with open(countFile, "rt") as fin:
        for line in fin:
          parts = line.rstrip().split('\t')
          if parts[1] == 'AA':
            totalAA += int(parts[2])

      logger.info("Total %d AA in dinucleotide file %s." % (totalAA, dinucleotideFile))

    logger.info("Processing dinucleotide file " + dinucleotideFile + " ...")
    
    dinuMap = finalMap[dinuName]

    count = 0
    tb = tabix.open(dinucleotideFile)
    for catItem in coordinates:
      count = count + 1
      if count % 10000 == 0:
        logger.info("%d / %d" % (count, len(coordinates)))

      if not catItem.reference_name in chromMap:
        continue
      
      #logger.info("Processing %s:%d-%d ..." %(catItem.reference_name, catItem.reference_start, catItem.reference_end))
      tbiter = tb.query(catItem.reference_name, catItem.reference_start, catItem.reference_end)
      records = [record for record in tbiter]
      for record in records:
        dinucleotide = record[3]
        diCount = int(record[4])
        if diCount == 0:
          diCount = 1

        start = int(record[1])
        index = catItem.getIndex(start)
        if index < 0:
          continue

        indexMap = dinuMap[index]
        
        if not dinucleotide in indexMap:
          indexMap[dinucleotide] = diCount
        else:
          indexMap[dinucleotide] = indexMap[dinucleotide] + diCount
    
  logger.info("Writing to %s ..." % outputFile)
  with open(outputFile, "wt") as fout:
    if normalizedByAA:
      fout.write("Sample\tPosition\tDinucleotide\tCount\tNormalizedCount\n")
    else:
      fout.write("Sample\tPosition\tDinucleotide\tCount\n")
    for dinuName in sorted( dinucleotideFileMap.keys() ):      
      for index in range(0, maxIndex):
        dinucleotideMap = finalMap[dinuName][index]
        for k in sorted(dinucleotideMap.keys()):
          count = dinucleotideMap[k]
          if normalizedByAA:
            normalizedCount = count * 1000000 / totalAA 
            fout.write("%s\t%s\t%s\t%d\t%lf\n" % (dinuName, index, k, count, normalizedCount))       
          else:
            fout.write("%s\t%s\t%s\t%d\n" % (dinuName, index, k, count))       

  rScript = os.path.join( os.path.dirname(__file__), "position.r")
  if not os.path.exists(rScript):
    raise Exception("Cannot find rscript %s" % rScript)
  
  cmd = "R --vanilla -f " + rScript + " --args " + os.path.abspath(outputFile) + " " + os.path.abspath(outputFile)
  runCmd(cmd, logger)

  logger.info("done.")
