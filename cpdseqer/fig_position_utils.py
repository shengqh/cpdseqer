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
from .common_utils import check_file_exists, get_reference_start, runCmd, readFileMap, checkFileMap, remove_chr, write_r_script

def readCoordinate(fin, addChr, delimit):
  result = []
  for line in fin:
    parts = line.rstrip().split(delimit)
    chrom = "chr" + parts[0] if addChr else parts[0] 
    result.append(CategoryItem(chrom, int(float(parts[1])), int(float(parts[2])), None))
  return(result)

def fig_position(logger, dinucleotideFileList, outputFile, coordinateFile, backgroundFile=None, useSpace=False, addChr=False, test=False):
  background_name = "__BACKGROUND__"

  logger.info("Reading dinucleotide file list ...")
  check_file_exists(dinucleotideFileList)
  dinucleotideFileMap = readFileMap(dinucleotideFileList)
  logger.info(str(dinucleotideFileMap))

  if (coordinateFile == 'hg38') or (coordinateFile == 'hg19'):
    coordinateFile = os.path.join(os.path.dirname(__file__), "data", "nucleosome_%s_interval.zip" % coordinateFile)
    useSpace = True
  check_file_exists(coordinateFile)

  rScript = os.path.join( os.path.dirname(__file__), "fig_position.r")
  if not os.path.exists(rScript):
    raise Exception("Cannot find rscript %s" % rScript)
  
  targetScript = write_r_script(outputFile, rScript, {'inputFile':outputFile})

  bNormalize = backgroundFile != None
  if bNormalize:
    check_file_exists(backgroundFile)
    dinucleotideFileMap[background_name] = backgroundFile

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
    countFile = dinucleotideFile.replace(".bed.bgz", ".count")

    check_file_exists(dinucleotideFile)
    check_file_exists(idxFile)
    check_file_exists(countFile)

    chromMap = {}
    with open(countFile, "rt") as fin:
      fin.readline()
      for line in fin:
        parts = line.split('\t')
        chromMap[parts[0]] = 1

    logger.info("Processing dinucleotide file %s : %s ..." % (dinuName, dinucleotideFile))
    
    dinuMap = finalMap[dinuName]

    count = 0
    tb = tabix.open(dinucleotideFile)
    for catItem in coordinates:
      count = count + 1
      if count % 10000 == 0:
        logger.info("%d / %d" % (count, len(coordinates)))

        if test and count % 100000 == 0:
          break

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
          indexMap[dinucleotide] = [1, diCount] # [siteCount, readCount]
        else:
          indexMap[dinucleotide][0] = indexMap[dinucleotide][0] + 1
          indexMap[dinucleotide][1] = indexMap[dinucleotide][1] + diCount
  
  if bNormalize:
    backgroundMap = finalMap.pop(background_name)

  logger.info("Writing to %s ..." % outputFile)
  with open(outputFile, "wt") as fout:
    if bNormalize:
      fout.write("Sample\tPosition\tDinucleotide\tSiteCount\tSiteNormalizedCount\tReadCount\tReadNormalizedCount\n")
    else:
      fout.write("Sample\tPosition\tDinucleotide\tSiteCount\tReadCount\n")
    for dinuName in sorted( finalMap.keys() ):      
      for index in range(0, maxIndex):
        dinucleotideMap = finalMap[dinuName][index]
        for k in sorted(dinucleotideMap.keys()):
          count = dinucleotideMap[k]
          if bNormalize:
            backgroundCount = backgroundMap[index][k] if k in backgroundMap[index] else [1, 1]
            fout.write("%s\t%s\t%s\t%d\t%lf\t%d\t%lf\n" % (dinuName, index, k, count[0], count[0] / backgroundCount[0], count[1], count[1] / backgroundCount[1]))       
          else:
            fout.write("%s\t%s\t%s\t%d\t%d\n" % (dinuName, index, k, count[0], count[1]))       

  options = {
    'inputFile':os.path.abspath(outputFile),
    'outfilePrefix':os.path.abspath(outputFile),
  }

  targetScript = write_r_script(outputFilePrefix, rScript, options)

  cmd = "R --vanilla -f " + targetScript
  runCmd(cmd, logger)

  logger.info("done.")
