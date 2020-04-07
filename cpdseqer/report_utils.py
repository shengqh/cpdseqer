import os
import os.path
import errno
import tabix
import sys
import shutil

from .CategoryItem import CategoryItem
from .common_utils import check_file_exists, get_reference_start, runCmd, readFileMap, checkFileMap, remove_chr

def report(logger, configFile, outputFilePrefix, block, dbVersion):
  check_file_exists(configFile)

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
  targetFolder = os.path.dirname(outputFilePrefix)
  targetConfigFile = os.path.join(targetFolder, "dinucleotide_file.list")
  with open(targetConfigFile, "wt") as fout:
    with open(configFile, "rt") as fin:
      for line in fin:
        parts = line.rstrip().split('\t')
        dinuFile = parts[0]
        dinuName = parts[1]
        dinuGroup = parts[2]
        targetDinuFile = os.path.join(targetFolder, "%s_block%d.txt" % (os.path.basename(dinuFile), block))
        fout.write("%s\t%s\t%s\n" % (targetDinuFile, dinuName, dinuGroup))

        if os.path.isfile(targetDinuFile):
          continue

        idxFile = dinuFile + ".tbi"

        check_file_exists(dinuFile)
        check_file_exists(idxFile)

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
          fdinu.write("%s\t%d\t%d\t%d\n" % (lastChrom, chromLengthMap[lastChrom] - 1, chromLengthMap[lastChrom], 0))

  targetRmd = outputFilePrefix + ".rmd"
  shutil.copyfile(rmdScript, targetRmd)

  cmd = "R -e \"setwd('%s');library(knitr);rmarkdown::render('%s',params=list(db='%s'));\"" % (os.path.dirname(targetRmd), os.path.basename(targetRmd), dbVersion)
  runCmd(cmd, logger)
