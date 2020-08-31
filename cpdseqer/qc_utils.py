import os
import os.path
import errno
import tabix
import sys
import shutil
import gzip

from .common_utils import check_file_exists, get_reference_start, runCmd, readFileMap, checkFileMap, remove_chr, write_rmd_script, ConfigItem, read_config_file
from .__version__ import __version__

def qc(logger, configFile, project_name, output_prefix, count_type, genome, size_factor_file):
  check_file_exists(configFile)

  items = read_config_file(configFile)
  logger.info("Config:" + str(items))

  if len(items) == 0:
    raise Exception("configFile %s is empty" % configFile)

  for item in items:
    check_file_exists(item.dinucleotide_file)
    check_file_exists(item.index_file)
    check_file_exists(item.count_file)
  
  if len(items) == 1:
    single_qc(logger, items[0], output_prefix, count_type, genome)
  else:
    multi_qc(logger, project_name, items, output_prefix, count_type, genome, size_factor_file)


def single_qc(logger, config_item, output_prefix, count_type, genome):
  options = {
    "sample": config_item.name,
    "bgz": config_item.dinucleotide_file,
    "cnt": config_item.count_file,
    "count_type": count_type,
    "gn": genome 
  }

  rScript = os.path.join( os.path.dirname(__file__), "qc.Rmd")
  targetScript = write_rmd_script(output_prefix, rScript, options)

  cmd = "R -e \"setwd('%s');library(knitr);rmarkdown::render('%s');\"" % (os.path.dirname(os.path.abspath(targetScript)), os.path.basename(targetScript))
  runCmd(cmd, logger)


def multi_qc(logger, project_name, config_items, output_prefix, count_type, genome, size_factor_file):
  targetConfigFile = output_prefix + ".config"
  with open(targetConfigFile, "wt") as fout:
    for item in config_items:
      fout.write("%s\t%s\t%s\n" % (item.name, item.count_file, item.dinucleotide_file))

  rScript = os.path.join( os.path.dirname(__file__), "multiQC.Rmd")
  
  options = {
    "smps":	os.path.basename(targetConfigFile),
    "dinuc": "DINUC4",
    "cntType": count_type,
    "exp": project_name,
    "gn": genome
  }

  if size_factor_file != None:
    options["libSizeNorm"] = size_factor_file

  targetScript = write_rmd_script(output_prefix, rScript, options)

  cmd = "R -e \"setwd('%s');library(knitr);rmarkdown::render('%s');\"" % (os.path.dirname(os.path.abspath(targetScript)), os.path.basename(targetScript))
  runCmd(cmd, logger)
