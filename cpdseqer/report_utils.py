import os
import os.path
import errno
import tabix
import sys
import shutil

from .CategoryItem import CategoryItem
from .common_utils import check_file_exists, get_reference_start, runCmd, readFileMap, checkFileMap, remove_chr
from .__version__ import __version__

def report(logger, sampleListFile, groupDefinitionFile, outputFilePrefix, block, dbVersion):
  check_file_exists(sampleListFile)
  check_file_exists(groupDefinitionFile)

  sampleMap = readFileMap(sampleListFile)
  groupMap = readFileMap(groupDefinitionFile)

  targetFolder = os.path.dirname(outputFilePrefix)

  checkFileMap(sampleMap)

  countFileMap = {sample:sampleMap[sample] + ".count" for sample in sampleMap.keys()}
  checkFileMap(countFileMap)

  sampleGroupMap = {}
  for sampleName in sampleMap.keys():
    if sampleName not in groupMap:
      raise Exception("Sample %s not in group definition file %s" % (sampleName, groupDefinitionFile))
    sampleGroupMap[sampleName] = [sampleMap[sampleName], groupMap[sampleName]]

  rmdScript = os.path.join( os.path.dirname(__file__), "report.rmd")
  if not os.path.exists(rmdScript):
    raise Exception("Cannot find rmdScript %s" % rmdScript)

  if os.path.isfile(dbVersion):
    chromInfo_file = dbVersion
  elif dbVersion == "hg19":
    chromInfo_file = os.path.join(os.path.dirname(__file__), "data", "chromInfo_hg19.txt")
  elif dbVersion == 'hg38':
    chromInfo_file = os.path.join(os.path.dirname(__file__), "data", "chromInfo_hg38.txt")
  else:
    raise Exception("I don't understand dbVersion: %s" % dbVersion)

  level_chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
  level_mut = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']

  logger.info("Reading chromosome length from: %s ..." % chromInfo_file)
  catItems = []# List<CategoryItem>
  chromLengthMap = {}
  with open(chromInfo_file, "rt") as fin:
    for line in fin:
      if line.startswith("#"):
        continue

      parts = line.rstrip().split('\t')
      chrom = parts[0]
      if not chrom in level_chr:
        continue

      chromLength = int(parts[1])
      chromLengthMap[chrom] = chromLength

      chromStart = 0
      while chromStart < chromLength:
        chromEnd = min(chromStart + block, chromLength)
        catItems.append(CategoryItem(chrom, chromStart, chromEnd, None))
        chromStart = chromStart + block

  logger.info("Preparing dinucleotide files ...")
  targetConfigFile = os.path.join(targetFolder, "dinucleotide_file.list")
  with open(targetConfigFile, "wt") as fout:
    fout.write("Group\tSample\tDinuFile\tCountFile\n")
    for dinuName in sorted(sampleGroupMap.keys()):
      dinuFile = sampleGroupMap[dinuName][0]
      dinuGroup = sampleGroupMap[dinuName][1]
      targetDinuFile = os.path.join(targetFolder, "%s_block%d.txt" % (os.path.basename(dinuFile), block))
      targetCountFile = os.path.join(targetFolder, "%s_count.txt" % os.path.basename(dinuFile))
      fout.write("%s\t%s\t%s\t%s\n" % (dinuGroup, dinuName, os.path.basename(targetDinuFile), os.path.basename(targetCountFile)))

      if os.path.isfile(targetDinuFile) and os.path.isfile(targetCountFile) :
        continue

      shutil.copyfile(countFileMap[dinuName], targetCountFile)

      idxFile = dinuFile + ".tbi"

      check_file_exists(dinuFile)
      check_file_exists(idxFile)

      dinuMap = {}
      logger.info("Preparing %s ..." % dinuFile)
      with open(targetDinuFile, "wt") as fdinu:
        tb = tabix.open(dinuFile)

        count = 0
        lastChrom = ""
        for catItem in catItems:
          count = count + 1
          if count % 10000 == 0:
            logger.info("%d / %d" % (count, len(catItems)))

          if catItem.reference_name != lastChrom:
            if lastChrom in chromLengthMap:
              fdinu.write("%s\t%d\t%d\t%d\n" % (lastChrom, chromLengthMap[lastChrom] - 1, chromLengthMap[lastChrom], 0))
            lastChrom = catItem.reference_name
          
          #logger.info("Processing %s:%d-%d ..." %(catItem.reference_name, catItem.reference_start, catItem.reference_end))
          tbiter = tb.query(catItem.reference_name, catItem.reference_start, catItem.reference_end)
          records = [record for record in tbiter]
          totalCount = 0
          for record in records:
            dinucleotide = record[3]
            diCount = int(record[4])
            if diCount == 0:
              diCount = 1
            totalCount += diCount

          if totalCount > 0:
            fdinu.write("%s\t%d\t%d\t%d\n" % (catItem.reference_name, catItem.reference_start, catItem.reference_end, totalCount))

  targetRmd = outputFilePrefix + ".rmd"
  shutil.copyfile(rmdScript, targetRmd)

  targetChromInfoFile =  os.path.join(targetFolder, "chromInfo.txt")
  shutil.copyfile(chromInfo_file, targetChromInfoFile)

  with open(os.path.join(targetFolder, "options.txt"), "wt") as fout:
    fout.write("Value\tName\n")
    fout.write("%s\tsoftVersion\n" % __version__)
    fout.write("%s\tdbVersion\n" % dbVersion)
    fout.write("%d\tblock\n" % block)

  cmd = "R -e \"setwd('%s');library(knitr);rmarkdown::render('%s');\"" % (os.path.dirname(targetRmd), os.path.basename(targetRmd))
  runCmd(cmd, logger)
