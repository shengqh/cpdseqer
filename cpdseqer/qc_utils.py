import os
import os.path
import errno
import tabix
import sys
import shutil

from .common_utils import check_file_exists, get_reference_start, runCmd, readFileMap, checkFileMap, remove_chr, write_rmd_script
from .__version__ import __version__

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

def qc(logger, sample_file_prefix, sample_name, output_file_prefix, count_type):
  targetFolder = os.path.dirname(output_file_prefix)

  dinuclotide_file = sample_file_prefix + ".bed.bgz"
  count_file = sample_file_prefix + ".count"
  check_file_exists(dinuclotide_file)
  check_file_exists(count_file)

  rScript = os.path.join( os.path.dirname(__file__), "qc.Rmd")
  
  options = {
    "sample": sample_name,
    "bgz": dinuclotide_file,
    "cnt": count_file,
    "count_type": count_type 
  }

  targetScript = write_rmd_script(output_file_prefix, rScript, options)

  cmd = "R -e \"setwd('%s');library(knitr);rmarkdown::render('%s');\"" % (os.path.dirname(os.path.abspath(targetScript)), os.path.basename(targetScript))
  runCmd(cmd, logger)
