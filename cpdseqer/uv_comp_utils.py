import os
import os.path
import errno
import tabix
import sys

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

def uv_comp_genome(logger, count_list_file, output_file_prefix, db, count_type):
  rScript = os.path.join( os.path.dirname(__file__), "stat_scenarios.Rmd")

  countFileMap = readFileMap(count_list_file)
  countFiles = ",".join([countFileMap[name] for name in sorted(countFileMap.keys())])

  options = {
    "scenario": "sce1",
    "smp": countFiles,
    "count_list_file": count_list_file,
    "gn": db,
    "cntType": count_type 
  }

  targetScript = write_rmd_script(output_file_prefix, rScript, options)

  cmd = "R -e \"setwd('%s');library(knitr);rmarkdown::render('%s');\"" % (os.path.dirname(os.path.abspath(targetScript)), os.path.basename(targetScript))
  runCmd(cmd, logger)


def uv_comp_genome_region(logger, count_list_file, output_file_prefix, db, count_type, coordinate_file):
  rScript = os.path.join( os.path.dirname(__file__), "stat_scenarios.Rmd")

  countFileMap = readFileMap(count_list_file)
  countFiles = ",".join([countFileMap[name] for name in sorted(countFileMap.keys())])

  options = {
    "scenario": "sce1",
    "smp": countFiles,
    "count_list_file": count_list_file,
    "gn": db,
    "cntType": count_type 
  }

  targetScript = write_rmd_script(output_file_prefix, rScript, options)

  cmd = "R -e \"setwd('%s');library(knitr);rmarkdown::render('%s');\"" % (os.path.dirname(os.path.abspath(targetScript)), os.path.basename(targetScript))
  runCmd(cmd, logger)
