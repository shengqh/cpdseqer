import os
import os.path
import errno
import tabix
import sys
import shutil
import gzip

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

def qc(logger, list_file, project_name, output_prefix, count_type):
  check_file_exists(list_file)

  files = []
  with open(list_file) as fh:
    for line in fh:
      parts = line.rstrip().split('\t')
      if len(parts) < 3:
        raise ArgumentError("list_file %s should contain at least three columns which indcates file name, count file and dinucleotide file without header")
      check_file_exists(parts[1])
      check_file_exists(parts[2])
      check_file_exists(parts[2] + ".tbi")
      files.append(parts)

  if len(files) == 0:
    raise ArgumentError("list_file doesn't is empty")

  if len(files) == 1:
    single_qc(logger, files[0][0], files[0][1], files[0][2], output_prefix, count_type)
  else:
    multi_qc(logger, project_name, list_file, output_prefix, count_type)


def single_qc(logger, sample_name, count_file, dinucleotide_file, output_prefix, count_type):
  targetFolder = os.path.dirname(output_prefix)

  options = {
    "sample": sample_name,
    "bgz": dinucleotide_file,
    "cnt": count_file,
    "count_type": count_type 
  }

  rScript = os.path.join( os.path.dirname(__file__), "qc.Rmd")
  targetScript = write_rmd_script(output_prefix, rScript, options)

  cmd = "R -e \"setwd('%s');library(knitr);rmarkdown::render('%s');\"" % (os.path.dirname(os.path.abspath(targetScript)), os.path.basename(targetScript))
  runCmd(cmd, logger)


def multi_qc(logger, project_name, list_file, output_prefix, count_type):
  targetFolder = os.path.dirname(output_prefix)

  rScript = os.path.join( os.path.dirname(__file__), "multiQC.Rmd")
  
  options = {
    "smps":	list_file,
    "dinuc": "DINUC4",
    "cntType": count_type,
    "exp": project_name
  }

  targetScript = write_rmd_script(output_prefix, rScript, options)

  cmd = "R -e \"setwd('%s');library(knitr);rmarkdown::render('%s');\"" % (os.path.dirname(os.path.abspath(targetScript)), os.path.basename(targetScript))
  runCmd(cmd, logger)
