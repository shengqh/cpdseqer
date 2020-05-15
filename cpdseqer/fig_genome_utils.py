import os
import os.path
import errno
import tabix
import sys
import shutil

from .CategoryItem import CategoryItem
from .common_utils import check_file_exists, get_reference_start, runCmd, readFileMap, checkFileMap, remove_chr, write_r_script
from .__version__ import __version__

class ConfigItem:
  def __init__(self, name, dinucleotide_file, is_case ):
    self.name = name
    self.dinucleotide_file = dinucleotide_file
    self.count_file = dinucleotide_file.replace(".bed.bgz", ".count")
    self.is_case = is_case
  
  def __repr__(self):
    return self.__str__()

  def __str__(self):
    return("[name=%s; dinucleotide_file=%s; count_file=%s; is_case=%s]" % (self.name, self.dinucleotide_file, self.count_file, self.is_case))

def read_config_file(config_file):
  result = []
  with open(config_file, "rt") as fin:
    headers = fin.readline().rstrip().split('\t')
    name_index = headers.index("Name")
    file_index = headers.index("File")
    case_index = headers.index("Case")
    for line in fin:
      parts = line.rstrip().split('\t')
      if len(parts) < len(headers):
        continue
      result.append(ConfigItem(parts[name_index], parts[file_index], parts[case_index]))
  return(result)

def fig_genome(logger, configFile, outputFilePrefix, block, dbVersion):
  check_file_exists(configFile)

  items = read_config_file(configFile)
  logger.info("Config:" + str(items))

  for item in items:
    check_file_exists(item.dinucleotide_file)
    check_file_exists(item.dinucleotide_file + ".tbi")
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
  level_mut = set(['CC', 'CT', 'TC', 'TT'])

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
    fout.write("Group\tSample\tDinuFile\tCountFile\n")
    for item in items:
      dinuName = item.name
      dinuGroup = item.is_case
      dinuFile = item.dinucleotide_file
      targetDinuFile = os.path.join(targetFolder, "%s_block%d.txt" % (os.path.basename(dinuFile), block))
      fout.write("%s\t%s\t%s\t%s\n" % (dinuGroup, dinuName, os.path.abspath(targetDinuFile), os.path.abspath(item.count_file)))

      #if os.path.isfile(targetDinuFile):
      #  continue

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
            if dinucleotide not in level_mut:
              continue

            diCount = int(record[4])
            if diCount == 0:
              diCount = 1
            totalReadCount += diCount
            totalSiteCount += 1

          if totalReadCount > 0:
            fdinu.write("%s\t%d\t%d\t%d\t%d\n" % (catItem.reference_name, catItem.reference_start, catItem.reference_end, totalReadCount, totalSiteCount))
        
  targetChromInfoFile =  os.path.join(targetFolder, "chromInfo.txt")
  shutil.copyfile(chromInfo_file, targetChromInfoFile)

  options = {
    'inputFile':targetConfigFile,
    'chromInfoFile':targetChromInfoFile,
    "block":block
  }

  targetScript = write_r_script(outputFilePrefix, rScript, options)

  cmd = "R --vanilla -f %s" % targetScript
  runCmd(cmd, logger)
