import os
import os.path
import errno
import tabix
import sys

from .DinucleotideItem import DinucleotideItem
from .CategoryItem import CategoryItem
from .common_utils import check_file_exists, get_reference_start, runCmd, readFileMap, checkFileMap, remove_chr

def statistic(logger, dinucleotideFileList, outputFile, coordinateFileList, category_index=-1, useSpace=False, addChr=False):
  check_file_exists(dinucleotideFileList)
  check_file_exists(coordinateFileList)

  dinucleotideFileMap = readFileMap(dinucleotideFileList)
  coordinateFileMap = readFileMap(coordinateFileList)
  
  logger.info("category_index=%d; useSpace=%s; addChr=%s" % (category_index, str(useSpace), str(addChr))) 

  coordinates = []
  delimit = ' ' if useSpace else '\t'
  for defCatName in coordinateFileMap.keys():
    coordinateFile = coordinateFileMap[defCatName]
    if (coordinateFile == 'hg38_promoter.bed') or (coordinateFile == 'hg38_tf.bed'):
      if not os.path.exists(coordinateFile):
        coordinateFile = os.path.join(os.path.dirname(__file__), "data", coordinateFile)

    check_file_exists(coordinateFile)
    
    logger.info("Reading category file " + coordinateFile + " ...")
      
    with open(coordinateFile, "rt") as fin:
      for line in fin:
        parts = line.rstrip().split(delimit)
        #print(parts)
        chrom = "chr" + parts[0] if addChr else parts[0] 
        catName = parts[category_index] if (category_index != -1 and category_index < len(parts)) else defCatName
        #print(catName)
        coordinates.append(CategoryItem(chrom, int(parts[1]), int(parts[2]), catName))
        #break
        
  catNames = sorted(list(set([ci.category for ci in coordinates])))

  finalMap = {dinuName:{catName:{} for catName in catNames} for dinuName in dinucleotideFileMap.keys()}
  
  for dinuName in dinucleotideFileMap.keys(): 
    dinucleotideFile = dinucleotideFileMap[dinuName]           
    idxFile = dinucleotideFile + ".tbi"

    check_file_exists(dinucleotideFile)
    check_file_exists(idxFile)

    logger.info("Processing dinucleotide file " + dinucleotideFile + " ...")
    
    count = 0
    tb = tabix.open(dinucleotideFile)
    for catItem in coordinates:
      count = count + 1
      if count % 10000 == 0:
        logger.info("%d / %d" % (count, len(coordinates)))
      
      #logger.info("Processing %s:%d-%d ..." %(catItem.reference_name, catItem.reference_start, catItem.reference_end))
      tbiter = tb.query(catItem.reference_name, catItem.reference_start, catItem.reference_end)
      records = [record for record in tbiter]
      for record in records:
        dinucleotide = record[3]
        count = int(record[4])
        if count == 0:
          count = 1

        if dinucleotide in catItem.dinucleotide_count_map.keys():
          catItem.dinucleotide_count_map[dinucleotide] = catItem.dinucleotide_count_map[dinucleotide] + count
        else:
          catItem.dinucleotide_count_map[dinucleotide] = count
    
    catDinucleotideMap = finalMap[dinuName]
    for ci in coordinates:
      dinucleotideMap = catDinucleotideMap[ci.category]
      for k in ci.dinucleotide_count_map.keys():
        if k in dinucleotideMap.keys():
          dinucleotideMap[k] = dinucleotideMap[k] + ci.dinucleotide_count_map[k]
        else:
          dinucleotideMap[k] = ci.dinucleotide_count_map[k]
  
  with open(outputFile, "wt") as fout:
    fout.write("Sample\tCategory\tDinucleotide\tCount\n")
    for dinuName in sorted( dinucleotideFileMap.keys() ):      
      for catName in catNames:
        dinucleotideMap = finalMap[dinuName][catName]
        for k in sorted(dinucleotideMap.keys()):
          fout.write("%s\t%s\t%s\t%d\n" % (dinuName, catName, k, dinucleotideMap[k]))    

  logger.info("done")   
