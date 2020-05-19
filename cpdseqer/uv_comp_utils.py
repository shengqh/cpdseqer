import os
import os.path
import errno
import tabix
import sys

from .common_utils import check_file_exists, check_data_file_exists, get_reference_start, runCmd, readFileMap, checkFileMap, remove_chr, write_rmd_script, md5sum
from .background_utils import genome_background, genome_background_region, calc_dinucleotide_distribution
from .count_utils import count
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

def calc_reg(logger, fasta_file):
  mut_levels = ['TT','TC','CC','CT']

  return []

def uv_comp_genome(logger, count_list_file, output_file_prefix, db, count_type):
  rScript = os.path.join( os.path.dirname(__file__), "stat_scenarios.Rmd")

  countFileMap = readFileMap(count_list_file)
  countFiles = ",".join([countFileMap[name] for name in sorted(countFileMap.keys())])

  options = {
    "scenario": "sce1",
    "smp": countFiles,
    "cntType": count_type 
  }

  targetFolder = os.path.dirname(os.path.abspath(output_file_prefix))

  if os.path.isfile(db):
    db_count_file = os.path.join(os.path.dirname(os.path.abspath(output_file_prefix)), os.path.basename(db) + ".count")
    if not os.path.isfile(db_count_file):
      logger.info(f"Build background table for {db}")
      genome_background(logger, db, db_count_file)
    reg = calc_dinucleotide_distribution(db_count_file)
    logger.info("background" + str(reg))
    options["reg"] = ",".join(reg)
  else:
    options["gb"] = db

  targetScript = write_rmd_script(output_file_prefix, rScript, options)

  cmd = "R -e \"setwd('%s');library(knitr);rmarkdown::render('%s');\"" % (targetFolder, os.path.basename(targetScript))
  runCmd(cmd, logger)

def uv_comp_genome_region(logger, dinu_list_file, output_file_prefix, db, count_type, coordinate_file, useSpace=False, addChr=False):

  #print(f"useSpace={useSpace}")

  check_file_exists(db)
  check_data_file_exists(coordinate_file)

  options = {
    "scenario": "sce1",
    "cntType": count_type 
  }

  coordinate_md5 = md5sum(coordinate_file)
  db_coord_count_file = os.path.join(os.path.dirname(os.path.abspath(output_file_prefix)), os.path.basename(db) + "." + coordinate_md5 + ".count")
  if not os.path.isfile(db_coord_count_file):
    logger.info(f"Build background table for {db}")
    genome_background_region(logger, db, db_coord_count_file, coordinate_file, useSpace, addChr)
  reg = calc_dinucleotide_distribution(db_coord_count_file)
  logger.info("background=" + str(reg))
  options["reg"] = ",".join(reg)

  sample_count_file = output_file_prefix + ".sample.count"
  count(logger, dinu_list_file, sample_count_file, coordinate_file, useSpace, addChr)
  options["smp"] = sample_count_file

  rScript = os.path.join( os.path.dirname(__file__), "stat_scenarios.Rmd")

  targetScript = write_rmd_script(output_file_prefix, rScript, options)

  cmd = "R -e \"setwd('%s');library(knitr);rmarkdown::render('%s');\"" % (os.path.dirname(os.path.abspath(targetScript)), os.path.basename(targetScript))
  runCmd(cmd, logger)
