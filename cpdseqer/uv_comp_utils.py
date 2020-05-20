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

def uv_comp_genome(logger, count_list_file, output_prefix, db, count_type):
  rScript = os.path.join( os.path.dirname(__file__), "stat_scenarios.Rmd")

  countFileMap = readFileMap(count_list_file)
  countFiles = ",".join([countFileMap[name] for name in sorted(countFileMap.keys())])

  options = {
    "scenario": "sce1",
    "smp": countFiles,
    "cntType": count_type 
  }

  targetFolder = os.path.dirname(os.path.abspath(output_prefix))

  if os.path.isfile(db):
    db_count_file = os.path.join(os.path.dirname(os.path.abspath(output_prefix)), os.path.basename(db) + ".count")
    if not os.path.isfile(db_count_file):
      logger.info(f"Build background table for {db}")
      genome_background(logger, db, db_count_file)
    reg = calc_dinucleotide_distribution(db_count_file)
    logger.info("background" + str(reg))
    options["reg"] = ",".join(reg)
  else:
    options["gb"] = db

  targetScript = write_rmd_script(output_prefix, rScript, options)

  cmd = "R -e \"setwd('%s');library(knitr);rmarkdown::render('%s');\"" % (targetFolder, os.path.basename(targetScript))
  runCmd(cmd, logger)

def calc_reg(logger, output_prefix, coordinate_file, db, useSpace, addChr):
  coordinate_md5 = md5sum(coordinate_file)
  db_coord_count_file = os.path.join(os.path.dirname(os.path.abspath(output_prefix)), os.path.basename(db) + "." + coordinate_md5 + ".count")
  if not os.path.isfile(db_coord_count_file):
    logger.info(f"Build background table for {db}")
    genome_background_region(logger, db, db_coord_count_file, coordinate_file, useSpace, addChr)
  reg = calc_dinucleotide_distribution(db_coord_count_file)
  return(reg)

def uv_comp_genome_region(logger, dinu_list_file, output_prefix, db, count_type, coordinate_file, useSpace=False, addChr=False):
  check_file_exists(db)
  coordinate_file = check_data_file_exists(coordinate_file)

  reg = calc_reg(logger, output_prefix, coordinate_file, db, useSpace, addChr)
  logger.info("background=" + str(reg))

  sample_count_file = output_prefix + ".sample.count"
  count(logger, dinu_list_file, sample_count_file, coordinate_file, useSpace, addChr)

  options = {
    "scenario": "sce1",
    "cntType": count_type ,
    "reg": ",".join(reg),
    "smp": sample_count_file
  }

  rScript = os.path.join( os.path.dirname(__file__), "stat_scenarios.Rmd")

  targetScript = write_rmd_script(output_prefix, rScript, options)

  cmd = "R -e \"setwd('%s');library(knitr);rmarkdown::render('%s');\"" % (os.path.dirname(os.path.abspath(targetScript)), os.path.basename(targetScript))
  runCmd(cmd, logger)

def uv_comp_regions(logger, dinu_list_file, output_prefix, db, count_type, coordinate_file1, coordinate_file2, useSpace=False, addChr=False):
  check_file_exists(db)
  coordinate_file1 = check_data_file_exists(coordinate_file1)
  coordinate_file2 = check_data_file_exists(coordinate_file2)

  options = {
    "scenario": "sce2",
    "cntType": count_type 
  }

  reg1 = calc_reg(logger, output_prefix, coordinate_file1, db, useSpace, addChr)
  logger.info(f"reg1={reg1}")

  reg2 = calc_reg(logger, output_prefix, coordinate_file2, db, useSpace, addChr)
  logger.info(f"reg2={reg2}")

  reg_count_file1 = output_prefix + ".reg1.count"
  count(logger, dinu_list_file, reg_count_file1, coordinate_file1, useSpace, addChr)

  reg_count_file2 = output_prefix + ".reg2.count"
  count(logger, dinu_list_file, reg_count_file2, coordinate_file2, useSpace, addChr)

  options = {
    "scenario": "sce2",
    "cntType": count_type ,
    "smp1": reg_count_file1,
    "smp2": reg_count_file2,
    "reg1": ",".join(reg1),
    "reg2": ",".join(reg2)
  }

  rScript = os.path.join( os.path.dirname(__file__), "stat_scenarios.Rmd")

  targetScript = write_rmd_script(output_prefix, rScript, options)

  cmd = "R -e \"setwd('%s');library(knitr);rmarkdown::render('%s');\"" % (os.path.dirname(os.path.abspath(targetScript)), os.path.basename(targetScript))
  runCmd(cmd, logger)

def uv_comp_groups(logger, count_list_file1, count_list_file2, output_prefix, count_type):
  rScript = os.path.join( os.path.dirname(__file__), "stat_scenarios.Rmd")

  countFileMap1 = readFileMap(count_list_file1)
  countFiles1 = ",".join([countFileMap1[name] for name in sorted(countFileMap1.keys())])

  countFileMap2 = readFileMap(count_list_file2)
  countFiles2 = ",".join([countFileMap2[name] for name in sorted(countFileMap2.keys())])

  options = {
    "scenario": "sce3",
    "smp1": countFiles1,
    "smp2": countFiles2,
    "cntType": count_type 
  }

  targetFolder = os.path.dirname(os.path.abspath(output_prefix))

  targetScript = write_rmd_script(output_prefix, rScript, options)

  cmd = "R -e \"setwd('%s');library(knitr);rmarkdown::render('%s');\"" % (targetFolder, os.path.basename(targetScript))
  runCmd(cmd, logger)


def uv_comp_groups_region(logger, dinu_list_file1, dinu_list_file2, output_prefix, count_type, coordinate_file, useSpace=False, addChr=False):
  rScript = os.path.join( os.path.dirname(__file__), "stat_scenarios.Rmd")

  coordinate_file = check_data_file_exists(coordinate_file)

  sample_count_file1 = output_prefix + ".sample1.count"
  count(logger, dinu_list_file1, sample_count_file1, coordinate_file, useSpace, addChr)

  sample_count_file2 = output_prefix + ".sample2.count"
  count(logger, dinu_list_file2, sample_count_file2, coordinate_file, useSpace, addChr)

  options = {
    "scenario": "sce3",
    "smp1": sample_count_file1,
    "smp2": sample_count_file2,
    "cntType": count_type 
  }

  targetFolder = os.path.dirname(os.path.abspath(output_prefix))

  targetScript = write_rmd_script(output_prefix, rScript, options)

  cmd = "R -e \"setwd('%s');library(knitr);rmarkdown::render('%s');\"" % (targetFolder, os.path.basename(targetScript))
  runCmd(cmd, logger)
