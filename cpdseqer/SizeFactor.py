import os
from collections import OrderedDict
from Bio import bgzf

from .common_utils import MUT_LEVELS, check_file_exists, write_r_script, read_config_file, runCmd

class AbstractSizeFactor:
  def __init__(self):
    self.suffix = ""

  def read_sites(self, logger, item):
    logger.info("Reading %s ..." % item.dinucleotide_file)
    curSites = {}
    line_count = 0
    with bgzf.BgzfReader(item.dinucleotide_file, "r") as fin:
      for line in fin:
        line_count += 1
        if line_count % 1000000 == 0:
          logger.info(line_count)
          
        parts = line.rstrip().split('\t')
        key = "%s_%s_%s" % (parts[0], parts[1], parts[3])
        curSites[key] = parts[4]
    return(curSites)

  def merge_sites(self, logger, total_sites, cur_sites, sample_name):
    raise Exception("No implementation")

  def calculate(self, logger, config_file, output_prefix):
    check_file_exists(config_file)

    items = read_config_file(config_file)
    logger.info("Config:" + str(items))

    output_folder = os.path.dirname(output_prefix)

    rScript = os.path.join( os.path.dirname(__file__), "size_factor.r")
    if not os.path.exists(rScript):
      raise Exception("Cannot find rScript %s" % rScript)

    bFirst = True
    sites = OrderedDict()

    for item in items:
      dinu_file = item.dinucleotide_file
      curSites = self.read_sites(logger, item)

      if bFirst:
        for site in curSites.keys():
          sites[site] = {item.name:curSites[site]}
        bFirst = False
        logger.info("Total %d sites readed from sample %s" % (len(sites), item.name))
        continue
      
      sites = self.merge_sites(logger, sites, curSites, item.name)

    target_count_file = output_prefix + self.suffix + ".count"
    with open(target_count_file, "wt") as fout:
      fout.write("Feature\t%s\n" % "\t".join([item.name for item in items]))
      for site in sites.keys():
        sample_dic = sites[site]
        fout.write("%s\t%s\n" % (site, "\t".join(sample_dic[item.name] if item.name in sample_dic else "0" for item in items)))

    output_file = output_prefix + self.suffix + ".sizefactor.txt"
    options = {
      'count_file':target_count_file,
      'output_file':output_file
    }

    targetScript = write_r_script(output_prefix + self.suffix, rScript, options)

    cmd = "R --vanilla -f %s" % targetScript
    runCmd(cmd, logger)


def intersect_features(logger, total_sites, cur_sites, sample_name):
  result = {}
  for site in cur_sites.keys():
    if site in total_sites:
      result[site] = total_sites[site]
      result[site][sample_name] = cur_sites[site]
  logger.info("Sites reduced to %d after intersected with sample %s" % (len(result), sample_name))
  return(result)

class SizeFactorSiteIntersection(AbstractSizeFactor):
  def __init__(self):
    self.suffix = ".site_intersect"

  def merge_sites(self, logger, total_sites, cur_sites, sample_name):
    return(intersect_features(logger, total_sites, cur_sites, sample_name))

def union_features(logger, total_sites, cur_sites, sample_name):
  result = total_sites
  for site in cur_sites.keys():
    result.setdefault(site, {})[sample_name] = cur_sites[site]
  logger.info("Sites increased to %d after unioned with sample %s" % (len(result), sample_name))
  return(result)

class SizeFactorSiteUnion(AbstractSizeFactor):
  def __init__(self):
    self.suffix = ".site_union"
    
  def merge_sites(self, logger, total_sites, cur_sites, sample_name):
    return (union_features(logger, total_sites, cur_sites, sample_name))

class SizeFactorChromDinucleotide(AbstractSizeFactor):
  def __init__(self):
    self.suffix = ".chrom_dinucleotide"

  def read_sites(self, logger, item):
    logger.info("Reading %s ..." % item.count_file)

    result = OrderedDict()
    with open(item.count_file, "rt") as fin:
      header = fin.readline()
      for line in fin:
        parts = line.rstrip().split('\t')

        #ignore the mutated dinu
        if parts[1] in MUT_LEVELS:
          continue

        #ignore the chromomsome contigs
        if len(parts[0]) > 4:
          continue

        key = parts[0] + "_" + parts[1]
        result[key] = parts[2]

    return(result)
    
  def merge_sites(self, logger, total_sites, cur_sites, sample_name):
    return (intersect_features(logger, total_sites, cur_sites, sample_name))
