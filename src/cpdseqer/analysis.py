import gzip
import pysam
import os
import tabix
from collections import OrderedDict
from Bio import SeqIO
from Bio import bgzf
from Bio.Seq import Seq

class DinucleotideItem:
  def __init__(self, reference_name, reference_start, reference_end, query_name, mapping_quality, strand, dinucleotide ):
    self.reference_name = reference_name
    self.reference_start = reference_start
    self.reference_end = reference_end
    self.query_name = query_name
    self.mapping_quality = mapping_quality
    self.strand = strand
    self.dinucleotide = dinucleotide;

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

def get_reference_start(elem):
    return elem.reference_start

def runCmd(cmd, logger):
  logger.info(cmd)
  os.system(cmd)

def readFileMap(fileName):
  result = {}
  with open(fileName) as fh:
    for line in fh:
      filepath, name = line.strip().split('\t', 1)
      result[name] = filepath.strip()
  return(result)
 
def demultiplex(logger, inputFile, outputFolder, configFile, args):
  logger.info("reading barcode file: " + configFile + " ...")
  sampleSeqMap = readFileMap(configFile)
  seqSampleMap = {sampleSeqMap[k]:k for k in sampleSeqMap.keys()}
  
  print(seqSampleMap)
  seqFileMap = {}
  barcodeLength = 0
  for seq in seqSampleMap.keys():
    barcodeLength = len(seq)
    seqFile = outputFolder + "/" + seqSampleMap[seq] + ".fastq.gz"
    seqFileMap[seq] = gzip.open(seqFile, 'wt')
  
  logger.info("reading input file: " + inputFile + " ...")
  count = 0
  with open(inputFile, "r") as fin:
    while(True):
      query = fin.readline()
      if not query:
        break
      
      seq = fin.readline().rstrip()
      skipline = fin.readline()
      score = fin.readline().rstrip()
      
      count = count + 1
      if count % 100000 == 0:
        logger.info("%d processed" % count)

      barcode = seq[0:barcodeLength]
      if barcode in seqFileMap:
        newseq = seq[barcodeLength:-1]
        newscore = score[barcodeLength:-1]
        fout = seqFileMap[barcode]
        fout.write(query)
        fout.write(newseq + '\n')
        fout.write(skipline)
        fout.write(newscore + '\n')
        
  for seq in seqFileMap.keys():
    seqFileMap[seq].close()
  
  logger.info("demultiplex done.")

def bam2dinucleotide(logger, bamFile, outputFile, genomeFastaFile):
  dinuItems = []
  count = 0
  with pysam.AlignmentFile(bamFile, "rb") as sf:
    for s in sf.fetch():
      count = count + 1
      if count % 100000 == 0:
        logger.info(count)
        
      if s.is_unmapped:
        continue
        
      if s.is_reverse:
        dinuItems.append(DinucleotideItem(s.reference_name, s.reference_end, s.reference_end + 2, s.query_name, s.mapping_quality, "-" if s.is_reverse else "+", "" ))
      else:
        dinuItems.append(DinucleotideItem(s.reference_name, s.reference_start - 2, s.reference_start, s.query_name, s.mapping_quality, "-" if s.is_reverse else "+", "" ))
  
  chrDinuMap = OrderedDict()
  for di in dinuItems:
    chrDinuMap.setdefault(di.reference_name, []).append(di)
  
  for values in chrDinuMap.values():
    values.sort(key=get_reference_start)
    
  with open(genomeFastaFile, "rt") as fin:  
    for record in SeqIO.parse(fin,'fasta'):
      id = record.id
      logger.info("Filling dinucleotide of " + id + " ...")

      if id in chrDinuMap.keys():
        seq = str(record.seq)
        seqlen = len(seq)
        chrDinuItems = chrDinuMap[id]
        for di in chrDinuItems:
          if di.reference_name == record.id:
            if di.reference_start >= 0 and di.reference_end <= seqlen:
              dinu = seq[di.reference_start:di.reference_end]
              if di.strand == "-":
                dinu = str(Seq(dinu).reverse_complement())
              di.dinucleotide = dinu
  
  logger.info("Writing dinucleotide to " + outputFile + " ...")
  with bgzf.BgzfWriter(outputFile, "wb") as fout:
    for chrom in chrDinuMap.keys():
      diList = chrDinuMap[chrom]
      for s in diList:
        if s.dinucleotide != "":
          fout.write("%s\t%d\t%d\t%s\t0\t%s\n" % (s.reference_name, s.reference_start, s.reference_end, s.dinucleotide, s.strand))
  
  runCmd("tabix -p bed %s " % outputFile, logger)
  logger.info("done.")
  
  #runCmd("bedtools flank -i %s -g %s -l 2 -r 0 -s > %s" % (bedFile, args.chromosome_size_file, dinuFile), logger)
  #runCmd("bedtools getfasta -fi %s -bed %s -s -fo %s" % (args.genome_seq_file, dinuFile, args.output), logger)
   
def remove_chr(chrom):
  if chrom.startswith("chr"):
    result = chrom[3:]
  else:
    result = chrom
  return(result)
  
def statistic(logger, dinucleotideFile, outputFile, coordinateFiles, coordinateFileNames = [], useSpace=False):
  coordinates = []
  hasName = len(coordinateFileNames) == len(coordinateFiles)
  delimit = ' ' if useSpace else '\t'
  for idx in range(0, len(coordinateFiles)):
    coordinateFile = coordinateFiles[idx]
    logger.info("Reading category file " + coordinateFile + " ...")
    if not hasName:
      defaultCat = os.path.splitext(os.path.basename(coordinateFile))[0]
    else:
      defaultCat = coordinateFileNames[idx]
      
    with open(coordinateFile, "rt") as fin:
      for line in fin:
        parts = line.rstrip().split(delimit)
        chrom = parts[0]
        catName = defaultCat if hasName else parts[3] if len(parts) >= 4 else defaultCat
        coordinates.append(CategoryItem("chr" + chrom, int(parts[1]), int(parts[2]), catName))
        
  logger.info("Processing dinucleotide file " + dinucleotideFile + " ...")

  count = 0
  tb = tabix.open(dinucleotideFile)
  for catItem in coordinates:
    count = count + 1
    if count % 10000 == 0:
      logger.info("%d / %d" % (count, len(coordinates)))
    
    #logger.info("Processing %s:%d-%d ..." %(catItem.reference_name, catItem.reference_start, catItem.reference_end))
    tbiter = tb.query(catItem.reference_name, catItem.reference_start, catItem.reference_end)
    records = [record for record in tbiter]
    for record in records:
      dinucleotide = record[3]
      if dinucleotide in catItem.dinucleotide_count_map.keys():
        catItem.dinucleotide_count_map[dinucleotide] = catItem.dinucleotide_count_map[dinucleotide] + 1
      else:
        catItem.dinucleotide_count_map[dinucleotide] = 1
  
  catDinucleotideMap = {}
  for ci in coordinates:
    if not ci.category in catDinucleotideMap.keys():
      catDinucleotideMap[ci.category] = {}
    dinucleotideMap = catDinucleotideMap[ci.category]
    for k in ci.dinucleotide_count_map.keys():
      if k in dinucleotideMap.keys():
        dinucleotideMap[k] = dinucleotideMap[k] + ci.dinucleotide_count_map[k]
      else:
        dinucleotideMap[k] = ci.dinucleotide_count_map[k]
  
  with open(outputFile, "wt") as fout:
    fout.write("Category\tDinucleotide\tCount\n")       
    for catName in sorted(catDinucleotideMap.keys()):
      dinucleotideMap = catDinucleotideMap[catName]
      for k in sorted(dinucleotideMap.keys()):
        if not 'N' in k:
          fout.write("%s\t%s\t%d\n" % (catName, k, dinucleotideMap[k]))       
  