import os
import os.path
import errno
import tabix
import sys

from .DinucleotideItem import DinucleotideItem
from .CategoryItem import CategoryItem
from .common_utils import check_file_exists, check_data_file_exists, get_reference_start, runCmd, readFileMap, checkFileMap, remove_chr, read_coordinate_file, read_chromosomes

def write_count_file(outputFile, finalMap):
  with open(outputFile, "wt") as fout:
    fout.write("Chromosome\tDinucleotide\tReadCount\tSiteCount\n")
    for catName in sorted( finalMap.keys() ):      
      dinuMap = finalMap[catName]
      for k in sorted(dinuMap.keys()):
        fout.write("%s\t%s\t%d\t%d\n" % (catName, k, dinuMap[k][0], dinuMap[k][1]))    

def count(logger, dinucleotideFileList, outputFile, coordinateFile, useSpace=False, addChr=False):
  check_file_exists(dinucleotideFileList)
  coordinateFile = check_data_file_exists(coordinateFile)

  dinucleotideFileMap = readFileMap(dinucleotideFileList)
  
  #logger.info("useSpace=%s; addChr=%s" % (str(useSpace), str(addChr))) 

  coordinates = []
  delimit = ' ' if useSpace else '\t'

  logger.info("Reading category file " + coordinateFile + " ...")
  coordinates = read_coordinate_file(coordinateFile, os.path.basename(coordinateFile), delimit, addChr)

  for coord in coordinates:
    coord.category = coord.reference_name
        
  catNames = sorted(list(set([ci.category for ci in coordinates])))

  finalMap = {catName:{} for catName in catNames}
  
  for dinuName in dinucleotideFileMap.keys(): 
    dinucleotideFile = dinucleotideFileMap[dinuName]           
    idxFile = dinucleotideFile + ".tbi"
    countFile = dinucleotideFile.replace(".bed.bgz", ".count")

    check_file_exists(dinucleotideFile)
    check_file_exists(idxFile)
    check_file_exists(countFile)

    chromosomes = read_chromosomes(countFile)

    logger.info("Processing dinucleotide file " + dinucleotideFile + " ...")
    
    catCount = 0
    tb = tabix.open(dinucleotideFile)
    for catItem in coordinates:
      catCount = catCount + 1
      if catCount % 10000 == 0:
        logger.info("%d / %d" % (catCount, len(coordinates)))

      if catItem.reference_name not in chromosomes:
        continue

      catMap = finalMap[catItem.reference_name]
      
      #logger.info("Processing %s:%d-%d ..." %(catItem.reference_name, catItem.reference_start, catItem.reference_end))
      tbiter = tb.query(catItem.reference_name, catItem.reference_start, catItem.reference_end)
      records = [record for record in tbiter]
      for record in records:
        dinu = record[3]
        count = int(record[4])
        if count == 0:
          count = 1

        if dinu in catMap.keys():
          catMap[dinu][0] += count
          catMap[dinu][1] += 1
        else:
          catMap[dinu] = [count, 1]
    
  write_count_file(outputFile, finalMap)

  logger.info("done")   
