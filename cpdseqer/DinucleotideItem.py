class DinucleotideItem:
  def __init__(self, reference_name, reference_start, reference_end, query_name, mapping_quality, strand, dinucleotide ):
    self.reference_name = reference_name
    self.reference_start = reference_start
    self.reference_end = reference_end
    self.query_name = query_name
    self.mapping_quality = mapping_quality
    self.strand = strand
    self.dinucleotide = dinucleotide
    self.count = 1
