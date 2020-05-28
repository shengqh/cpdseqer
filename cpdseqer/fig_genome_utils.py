import os
import os.path
import errno
import tabix
import sys
import shutil

from .CategoryItem import CategoryItem
from .common_utils import MUT_LEVELS, check_file_exists, get_reference_start, runCmd, readFileMap, checkFileMap, remove_chr, write_r_script, ConfigItem, read_config_file
from .__version__ import __version__

def fig_genome(logger, configFile, outputFilePrefix, block, dbVersion, raw_count):
  check_file_exists(configFile)

  items = read_config_file(configFile)
  logger.info("Config:" + str(items))

  for item in items:
    check_file_exists(item.dinucleotide_file)
    check_file_exists(item.index_file)
    check_file_exists(item.count_file)

  targetFolder = os.path.dirname(outputFilePrefix)

  rScript = os.path.join( os.path.dirname(__file__), "fig_genome.r")
  if not os.path.exists(rScript):
    raise Exception("Cannot find rScript %s" % rScript)

  if os.path.isfile(dbVersion):
    chromInfo_file = dbVersion
  elif dbVersion == "hg19":
    chromInfo_file = os.path.join(os.path.dirname(__file__), "data", "chromInfo_hg19.txt")
  elif dbVersion == 'hg38':
    chromInfo_file = os.path.join(os.path.dirname(__file__), "data", "chromInfo_hg38.txt")
  else:
    raise Exception("I don't understand dbVersion: %s" % dbVersion)

  level_chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
  #level_mut = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
  level_mut = set(MUT_LEVELS)

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

  logger.info("Processing dinucleotide files ...")
  targetConfigFile = outputFilePrefix + ".config.txt"
  with open(targetConfigFile, "wt") as fout:
    fout.write("Sample\tDinuFile\tCountFile\tTotalRead\tTotalSite\n")
    for item in items:
      dinuName = item.name
      dinuFile = item.dinucleotide_file
      targetDinuFile = os.path.join(targetFolder, "%s_block%d.txt" % (os.path.basename(dinuFile), block))

      #if os.path.isfile(targetDinuFile):
      #  continue

      totalFileReadCount = 0
      totalFileSiteCount = 0
      dinuMap = {}
      logger.info("Processing %s ..." % dinuFile)
      with open(targetDinuFile, "wt") as fdinu:
        fdinu.write("Chrom\tStart\tEnd\tReadCount\tSiteCount\n")
        tb = tabix.open(dinuFile)

        count = 0
        for catItem in catItems:
          count = count + 1
          if count % 10000 == 0:
            logger.info("%d / %d" % (count, len(catItems)))

          #logger.info("Processing %s:%d-%d ..." %(catItem.reference_name, catItem.reference_start, catItem.reference_end))
          tbiter = tb.query(catItem.reference_name, catItem.reference_start, catItem.reference_end)
          records = [record for record in tbiter]
          totalReadCount = 0
          totalSiteCount = 0
          for record in records:
            dinucleotide = record[3]

            diCount = int(record[4])
            if diCount == 0:
              diCount = 1

            totalFileReadCount += diCount
            totalFileSiteCount += 1

            if dinucleotide not in level_mut:
              continue

            totalReadCount += diCount
            totalSiteCount += 1

          if totalReadCount > 0:
            fdinu.write("%s\t%d\t%d\t%d\t%d\n" % (catItem.reference_name, catItem.reference_start, catItem.reference_end, totalReadCount, totalSiteCount))

      fout.write("%s\t%s\t%s\t%d\t%d\n" % (dinuName, os.path.abspath(targetDinuFile), os.path.abspath(item.count_file), totalFileReadCount, totalFileSiteCount))
        
  targetChromInfoFile =  os.path.join(targetFolder, "chromInfo.txt")
  shutil.copyfile(chromInfo_file, targetChromInfoFile)

  options = {
    'inputFile':targetConfigFile,
    'chromInfoFile':targetChromInfoFile,
    "block":block,
    'useRawCount':"1" if raw_count else "0"
  }

  targetScript = write_r_script(outputFilePrefix, rScript, options)

  cmd = "R --vanilla -f %s" % targetScript
  runCmd(cmd, logger)
