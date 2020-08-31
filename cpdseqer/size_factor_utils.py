from .SizeFactor import SizeFactorSiteIntersection, SizeFactorSiteUnion, SizeFactorChromDinucleotide

def size_factor(logger, config_file, output_prefix, calc_type):
  #if calc_type == "site_intersect":
  #  sf_calc = SizeFactorSiteIntersection()
  #elif 
  if calc_type == "site_union":
    sf_calc = SizeFactorSiteUnion()
  elif calc_type == "chrom_dinucleotide":
    sf_calc = SizeFactorChromDinucleotide()
  else: 
    raise Exception("Not implemented for calc_type %s" % calc_type)
  
  sf_calc.calculate(logger, config_file, output_prefix)
