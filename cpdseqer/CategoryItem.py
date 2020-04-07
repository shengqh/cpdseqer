class CategoryItem:
  def __init__(self, reference_name, reference_start, reference_end, category ):
    self.reference_name = reference_name
    self.reference_start = reference_start
    self.reference_end = reference_end
    self.category = category
    self.dinucleotide_count_map = {}
    
  def overlap (self, position):
    if position < self.reference_start:
      return(1)
    if position <= self.reference_end:
      return(0)
    return(-1)

  def getIndex(self, position):
    return(position - self.reference_start)

  def getLength(self):
    return(self.reference_end - self.reference_start)
