import gzip
import zipfile
import pysam
import os
import os.path
import errno
import tabix
import sys
import shutil
from io import TextIOWrapper

from .CategoryItem import CategoryItem
from .common_utils import check_file_exists, get_reference_start, runCmd, readFileMap, checkFileMap, remove_chr


def readCoordinate(fin, addChr, delimit):
  result = []
  for line in fin:
    parts = line.rstrip().split(delimit)
    chrom = "chr" + parts[0] if addChr else parts[0] 
    result.append(CategoryItem(chrom, int(float(parts[1])), int(float(parts[2])), None))
  return(result)

def position(logger, dinucleotideFileList, outputFile, coordinateFile, normalizedFiles=None, useSpace=False, addChr=False):
  check_file_exists(dinucleotideFileList)
  dinucleotideFileMap = readFileMap(dinucleotideFileList)

  if (coordinateFile == 'hg38') or (coordinateFile == 'hg19'):
    coordinateFile = os.path.join(os.path.dirname(__file__), "data", "nucleosome_%s_interval.zip" % coordinateFile)
  check_file_exists(coordinateFile)

  bNormalize = normalizedFiles != None
  if bNormalize:
    if normalizedFiles == "auto":
      normFiles = [os.path.join(os.path.dirname(__file__), "data", 'Naked.1.count'), os.path.join(os.path.dirname(__file__), "data", 'Naked.2.count')]
    else:
      normFiles = normalizedFiles.split(',')

    for normFile in normFiles:
      check_file_exists(normFile)

    backgroundMap = {}
    for normFile in normFiles:
      with open(normFile, "rt") as fin:
        fin.readline()
        for line in fin:
          parts = line.rstrip().split('\t')
          if parts[1] in backgroundMap:
            backgroundMap[parts[1]] += int(parts[2])
          else:
            backgroundMap[parts[1]] = int(parts[2])
    
    meanCount = sum(backgroundMap.values()) / len(backgroundMap)
    backgroundFactor = { key:(backgroundMap[key] / meanCount) for key in backgroundMap.keys() }
    with open(outputFile + ".norm.txt", "wt") as fout:
      fout.write("Dinucleotide\tCount\tFactor\n")
      for dinu in sorted(backgroundMap.keys()):
        fout.write("%s\t%d\t%.lf" % (dinu, backgroundMap[dinu], backgroundFactor[dinu]))

  logger.info("Reading coordinate file " + coordinateFile + " ...")
      
  delimit = ' ' if useSpace else '\t'
  if coordinateFile.endswith(".zip"):
    internalFile = os.path.basename(coordinateFile).replace(".zip", ".txt")
    with zipfile.ZipFile(coordinateFile) as z:
      with z.open(internalFile) as fz:
        with TextIOWrapper(fz) as fin:
          coordinates = readCoordinate(fin, addChr, delimit)
  elif coordinateFile.endswith(".gz"):
    with gzip.open(coordinateFile, "rt") as fin:
      coordinates = readCoordinate(fin, addChr, delimit)
  else:
    with open(coordinate, "rt") as fin:
      coordinates = readCoordinate(fin, addChr, delimit)

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
    if bNormalize:
      fout.write("Sample\tPosition\tDinucleotide\tCount\tNormalizedCount\n")
    else:
      fout.write("Sample\tPosition\tDinucleotide\tCount\n")
    for dinuName in sorted( dinucleotideFileMap.keys() ):      
      for index in range(0, maxIndex):
        dinucleotideMap = finalMap[dinuName][index]
        for k in sorted(dinucleotideMap.keys()):
          count = dinucleotideMap[k]
          if bNormalize:
            normalizedCount = count / backgroundFactor[k] 
            fout.write("%s\t%s\t%s\t%d\t%lf\n" % (dinuName, index, k, count, normalizedCount))       
          else:
            fout.write("%s\t%s\t%s\t%d\n" % (dinuName, index, k, count))       

  rScript = os.path.join( os.path.dirname(__file__), "position.r")
  if not os.path.exists(rScript):
    raise Exception("Cannot find rscript %s" % rScript)
  
  cmd = "R --vanilla -f " + rScript + " --args " + os.path.abspath(outputFile) + " " + os.path.abspath(outputFile)
  runCmd(cmd, logger)

  logger.info("done.")
