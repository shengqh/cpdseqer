import os
import os.path
import logging
import errno
import shutil
import hashlib

MUT_LEVELS=['TT','TC','CC','CT']
DINU_LEVELS=['AA', 'AC', 'AT', 'AG', 'CC', 'CA', 'CT', 'CG', 'GG', 'GA', 'GC', 'GT','TT', 'TA', 'TC', 'TG']

from .CategoryItem import CategoryItem

def check_file_exists(file):
  if not os.path.exists(file):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)

def check_data_file_exists(file):
  if os.path.exists(file):
    return(file)
  
  possibleFile = os.path.join(os.path.dirname(__file__), "data", file)
  if os.path.exists(possibleFile):
    return possibleFile
  
  raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)

def get_reference_start(elem):
    return elem.reference_start

def runCmd(cmd, logger):
  logger.info(cmd)
  os.system(cmd)

def readFileMap(fileName):
  check_file_exists(fileName)

  result = {}
  with open(fileName) as fh:
    for line in fh:
      filepath, name = line.strip().split('\t', 1)
      result[name] = filepath.strip()
  return(result)

def checkFileMap(fileMap):
  for sname in fileMap.keys():
    sfile = fileMap[sname]
    check_file_exists(sfile)

def remove_chr(chrom):
  if chrom.startswith("chr"):
    result = chrom[3:]
  else:
    result = chrom
  return(result)

def initialize_logger(logfile, args):
  logger = logging.getLogger('cpdseqer')
  loglevel = logging.INFO
  logger.setLevel(loglevel)

  formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')    
 
  # create console handler and set level to info
  handler = logging.StreamHandler()
  handler.setLevel(loglevel)
  handler.setFormatter(formatter)
  logger.addHandler(handler)
 
  # create error file handler and set level to error
  handler = logging.FileHandler(logfile, "w")
  handler.setLevel(loglevel)
  handler.setFormatter(formatter)
  logger.addHandler(handler)
 
  return(logger)

def read_coordinate_file(fileName, defCatName, delimit='\t', addChr=False, categoryIndex=-1):
  #print("delimit=tab" if delimit=='\t' else "delimit=space")
  result = []
  with open(fileName, "rt") as fin:
    for line in fin:
      parts = line.rstrip().split(delimit)
      #print(parts)
      chrom = "chr" + parts[0] if addChr else parts[0] 
      catName = parts[categoryIndex] if (categoryIndex != -1 and categoryIndex < len(parts)) else defCatName
      #print(catName)
      result.append(CategoryItem(chrom, int(float(parts[1])), int(float(parts[2])), catName))
      
  return(result)

def write_r_script(outfilePrefix, rScript, optionMap={}):
  targetScript = outfilePrefix + ".r"
  optionMap["outfilePrefix"] = outfilePrefix
  with open(targetScript, "wt") as fout:
    for key in optionMap.keys():
      fout.write("%s='%s'\n" % (key, optionMap[key]))

    fout.write("setwd('%s')\n" % os.path.dirname(os.path.abspath(targetScript)))
    
    with open(rScript, "rt") as fin:
      bFirstSetwd = True
      for line in fin:
        if line.startswith("setwd") and bFirstSetwd:
          bFirstSetwd = False
          continue

        bInOption = False
        for key in optionMap.keys():
          if line.startswith(key + "="):
            optionMap.pop(key)
            bInOption = True
            break

        if not bInOption:
          fout.write(line)

  return(targetScript)

def write_rmd_script(outfilePrefix, rmdScript, optionMap={}, copyRFunction=True, optionsToIndividualFile=True):
  if not os.path.exists(rmdScript):
    raise Exception("Cannot find file %s" % rmdScript)

  if copyRFunction:
    rFunScript = os.path.join( os.path.dirname(__file__), "Rfunctions.R")
    if not os.path.exists(rFunScript):
      raise Exception("Cannot find rScript %s" % rFunScript)

    targetFolder = os.path.dirname(os.path.abspath(outfilePrefix))
    targetRFunScript =  os.path.join(targetFolder, "Rfunctions.R")
    shutil.copyfile(rFunScript, targetRFunScript)

  if optionsToIndividualFile:
    optionToIndividualFileName = outfilePrefix + ".options"
    with open(optionToIndividualFileName, "wt") as fout:
      for key in sorted(optionMap.keys()):
        fout.write("%s\t%s\n" % (key, optionMap[key]))
    optionMap = {"option_file": os.path.basename(optionToIndividualFileName)}

  targetScript = outfilePrefix + ".rmd"
  with open(targetScript, "wt") as fout:
    with open(rmdScript, "rt") as fin:
      for line in fin:
        if line.startswith("```"):
          fout.write(line)
          for key in optionMap.keys():
            fout.write("%s='%s'\n" % (key, optionMap[key]))
          fout.write("\n")
          break
        else:
          fout.write(line)

      for line in fin:
        bInOption = False
        for key in optionMap.keys():
          if line.startswith(key + "=") or line.startswith(key + "<-"):
            optionMap.pop(key)
            bInOption = True
            break

        if not bInOption:
          fout.write(line)

  return(targetScript)

def md5sum(filename, blocksize=65536):
    hash = hashlib.md5()
    with open(filename, "rb") as f:
        for block in iter(lambda: f.read(blocksize), b""):
            hash.update(block)
    return hash.hexdigest()

def read_chromosomes(countFile):
    chromMap = {}
    with open(countFile, "rt") as fin:
      fin.readline()
      for line in fin:
        parts = line.split('\t')
        chromMap[parts[0]] = 1
    return (sorted(list(chromMap.keys())))
    
def get_count_file(dinucleotide_file):
  return(dinucleotide_file.replace(".bed.gz", ".count"))
