import gzip
import pysam
import os
import os.path
import errno
import tabix
import sys
import shutil
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
    self.dinucleotide = dinucleotide
    self.count = 1

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


def check_file_exists(file):
  if not os.path.exists(file):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)

def get_reference_start(elem):
    return elem.reference_start

def runCmd(cmd, logger):
  logger.info(cmd)
  os.system(cmd)

def readFileMap(fileName):
  check_file_exists(fileName)

  result = {}
  with open(fileName) as fh:
    for line in fh:
      filepath, name = line.strip().split('\t', 1)
      result[name] = filepath.strip()
  return(result)

def checkFileMap(fileMap):
  for sname in fileMap.keys():
    sfile = fileMap[sname]
    check_file_exists(sfile)
 
def demultiplex(logger, inputFile, outputFolder, configFile, args):
  check_file_exists(inputFile)

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
  
  if inputFile.endswith(".gz"):
    fin = gzip.open(inputFile, 'rt')
  else:
    fin = open(inputFile, "rt")

  with fin:
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

def bam2dinucleotide(logger, bamFile, outputFile, genomeFastaFile, mappingQuality=20):
  check_file_exists(bamFile)
  check_file_exists(genomeFastaFile)

  logger.info("reading bam file %s ..." % bamFile )
  dinuItems = []
  count = 0
  with pysam.AlignmentFile(bamFile, "rb") as sf:
    for s in sf.fetch():
      count = count + 1
      if count % 100000 == 0:
        logger.info(count)
        
      if s.is_unmapped:
        continue

      if s.mapping_quality < mappingQuality:
        continue
        
      if s.is_reverse:
        dinuItems.append(DinucleotideItem(s.reference_name, s.reference_end, s.reference_end + 2, s.query_name, s.mapping_quality, "-", "" ))
      else:
        dinuItems.append(DinucleotideItem(s.reference_name, s.reference_start - 2, s.reference_start, s.query_name, s.mapping_quality, "+", "" ))
  
  chrDinuMap = OrderedDict()
  for di in dinuItems:
    chrDinuMap.setdefault(di.reference_name, []).append(di)
  
  for chr in chrDinuMap.keys():
    values = chrDinuMap[chr]
    logger.info("sort %d dinucleotides of chromosome %s..." % (len(values), chr ))
    values.sort(key=get_reference_start)
    logger.info("combine %d dinucleotides of chromosome %s..." % (len(values), chr ) )
    idx = len(values) - 1
    deleteList = set()
    while(idx > 0):
      curDinu = values[idx]
      prev = idx -1
      while(prev >= 0):
        prevDinu = values[prev]
        if curDinu.reference_start != prevDinu.reference_start:
          break
        if curDinu.strand == prevDinu.strand:
          prevDinu.count = prevDinu.count + curDinu.count
          deleteList.add(idx)
          break
        prev = prev - 1
      idx = idx -1
    chrDinuMap[chr] = [i for j, i in enumerate(values) if j not in deleteList]
    logger.info("after combine, there is %d dinucleotides of chromosome %s..." % (len(chrDinuMap[chr]), chr ) )

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
              if di.strand == "+":
                dinu = str(Seq(dinu).reverse_complement())
              di.dinucleotide = dinu
  
  countMap = OrderedDict()
  logger.info("Writing dinucleotide to " + outputFile + " ...")
  with bgzf.BgzfWriter(outputFile, "wb") as fout:
    for chrom in chrDinuMap.keys():
      diList = chrDinuMap[chrom]
      chromMap = {}
      for s in diList:
        if (s.dinucleotide != "") and (not 'N' in s.dinucleotide):
          chromMap[s.dinucleotide] = chromMap.setdefault(s.dinucleotide, 0) + s.count
          fout.write("%s\t%d\t%d\t%s\t%d\t%s\n" % (s.reference_name, s.reference_start, s.reference_end, s.dinucleotide, s.count, s.strand))
      countMap[chrom] = chromMap
  
  with open(outputFile + ".count", "wt") as fout:
    fout.write("Chromosome\tDinucleotide\tCount\n")
    for chrom in countMap.keys():
      chromMap = countMap[chrom]
      for dinucleotide in chromMap.keys():
        fout.write("%s\t%s\t%d\n" % (chrom, dinucleotide, chromMap[dinucleotide]))

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
  
def statistic(logger, dinucleotideFileList, outputFile, coordinateFileList, category_index=-1, useSpace=False, addChr=False):
  check_file_exists(dinucleotideFileList)
  check_file_exists(coordinateFileList)

  dinucleotideFileMap = readFileMap(dinucleotideFileList)
  coordinateFileMap = readFileMap(coordinateFileList)
  
  logger.info("category_index=%d; useSpace=%s; addChr=%s" % (category_index, str(useSpace), str(addChr))) 

  coordinates = []
  delimit = ' ' if useSpace else '\t'
  for defCatName in coordinateFileMap.keys():
    coordinateFile = coordinateFileMap[defCatName]
    if (coordinateFile == 'hg38_promoter.bed') or (coordinateFile == 'hg38_tf.bed'):
      if not os.path.exists(coordinateFile):
        coordinateFile = os.path.join(os.path.dirname(__file__), "data", coordinateFile)

    check_file_exists(coordinateFile)
    
    logger.info("Reading category file " + coordinateFile + " ...")
      
    with open(coordinateFile, "rt") as fin:
      for line in fin:
        parts = line.rstrip().split(delimit)
        #print(parts)
        chrom = "chr" + parts[0] if addChr else parts[0] 
        catName = parts[category_index] if (category_index != -1 and category_index < len(parts)) else defCatName
        #print(catName)
        coordinates.append(CategoryItem(chrom, int(parts[1]), int(parts[2]), catName))
        #break
        
  catNames = sorted(list(set([ci.category for ci in coordinates])))

  finalMap = {dinuName:{catName:{} for catName in catNames} for dinuName in dinucleotideFileMap.keys()}
  
  for dinuName in dinucleotideFileMap.keys(): 
    dinucleotideFile = dinucleotideFileMap[dinuName]           
    idxFile = dinucleotideFile + ".tbi"

    check_file_exists(dinucleotideFile)
    check_file_exists(idxFile)

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
        count = int(record[4])
        if count == 0:
          count = 1

        if dinucleotide in catItem.dinucleotide_count_map.keys():
          catItem.dinucleotide_count_map[dinucleotide] = catItem.dinucleotide_count_map[dinucleotide] + count
        else:
          catItem.dinucleotide_count_map[dinucleotide] = count
    
    catDinucleotideMap = finalMap[dinuName]
    for ci in coordinates:
      dinucleotideMap = catDinucleotideMap[ci.category]
      for k in ci.dinucleotide_count_map.keys():
        if k in dinucleotideMap.keys():
          dinucleotideMap[k] = dinucleotideMap[k] + ci.dinucleotide_count_map[k]
        else:
          dinucleotideMap[k] = ci.dinucleotide_count_map[k]
  
  with open(outputFile, "wt") as fout:
    fout.write("Sample\tCategory\tDinucleotide\tCount\n")
    for dinuName in sorted( dinucleotideFileMap.keys() ):      
      for catName in catNames:
        dinucleotideMap = finalMap[dinuName][catName]
        for k in sorted(dinucleotideMap.keys()):
          fout.write("%s\t%s\t%s\t%d\n" % (dinuName, catName, k, dinucleotideMap[k]))       

def position(logger, dinucleotideFileList, outputFile, coordinateFile, normalizedByAA=False, useSpace=False, addChr=False):
  check_file_exists(dinucleotideFileList)
  check_file_exists(coordinateFile)

  dinucleotideFileMap = readFileMap(dinucleotideFileList)

  coordinates = []
  delimit = ' ' if useSpace else '\t'
  if (coordinateFile == 'hg38_promoter.bed') or (coordinateFile == 'hg38_tf.bed'):
    if not os.path.exists(coordinateFile):
      coordinateFile = os.path.join(os.path.dirname(__file__), "data", coordinateFile)

  logger.info("Reading coordinate file " + coordinateFile + " ...")
      
  with open(coordinateFile, "rt") as fin:
    for line in fin:
      parts = line.rstrip().split(delimit)
      chrom = "chr" + parts[0] if addChr else parts[0] 
      coordinates.append(CategoryItem(chrom, int(float(parts[1])), int(float(parts[2])), None))

  maxIndex = max([cor.getLength() for cor in coordinates])
        
  finalMap = {dinuName:{index:{} for index in range(0, maxIndex)} for dinuName in dinucleotideFileMap.keys()}
  
  for dinuName in dinucleotideFileMap.keys(): 
    dinucleotideFile = dinucleotideFileMap[dinuName]           
    idxFile = dinucleotideFile + ".tbi"
    countFile = dinucleotideFile + ".count"

    check_file_exists(dinucleotideFile)
    check_file_exists(idxFile)
    check_file_exists(countFile)

    chromMap = {}
    with open(countFile, "rt") as fin:
      fin.readline()
      for line in fin:
        parts = line.split('\t')
        chromMap[parts[0]] = 1

    totalAA = 0
    if normalizedByAA:
      countFile = dinucleotideFile + ".count"
      check_file_exists(countFile)

      with open(countFile, "rt") as fin:
        for line in fin:
          parts = line.rstrip().split('\t')
          if parts[1] == 'AA':
            totalAA += int(parts[2])

      logger.info("Total %d AA in dinucleotide file %s." % (totalAA, dinucleotideFile))

    logger.info("Processing dinucleotide file " + dinucleotideFile + " ...")
    
    dinuMap = finalMap[dinuName]

    count = 0
    tb = tabix.open(dinucleotideFile)
    for catItem in coordinates:
      count = count + 1
      if count % 10000 == 0:
        logger.info("%d / %d" % (count, len(coordinates)))

      if not catItem.reference_name in chromMap:
        continue
      
      #logger.info("Processing %s:%d-%d ..." %(catItem.reference_name, catItem.reference_start, catItem.reference_end))
      tbiter = tb.query(catItem.reference_name, catItem.reference_start, catItem.reference_end)
      records = [record for record in tbiter]
      for record in records:
        dinucleotide = record[3]
        diCount = int(record[4])
        if diCount == 0:
          diCount = 1

        start = int(record[1])
        index = catItem.getIndex(start)
        if index < 0:
          continue

        indexMap = dinuMap[index]
        
        if not dinucleotide in indexMap:
          indexMap[dinucleotide] = diCount
        else:
          indexMap[dinucleotide] = indexMap[dinucleotide] + diCount
    
  logger.info("Writing to %s ..." % outputFile)
  with open(outputFile, "wt") as fout:
    if normalizedByAA:
      fout.write("Sample\tPosition\tDinucleotide\tCount\tNormalizedCount\n")
    else:
      fout.write("Sample\tPosition\tDinucleotide\tCount\n")
    for dinuName in sorted( dinucleotideFileMap.keys() ):      
      for index in range(0, maxIndex):
        dinucleotideMap = finalMap[dinuName][index]
        for k in sorted(dinucleotideMap.keys()):
          count = dinucleotideMap[k]
          if normalizedByAA:
            normalizedCount = count * 1000000 / totalAA 
            fout.write("%s\t%s\t%s\t%d\t%lf\n" % (dinuName, index, k, count, normalizedCount))       
          else:
            fout.write("%s\t%s\t%s\t%d\n" % (dinuName, index, k, count))       

  rScript = os.path.join( os.path.dirname(__file__), "position.r")
  if not os.path.exists(rScript):
    raise Exception("Cannot find rscript %s" % rScript)
  
  cmd = "R --vanilla -f " + rScript + " --args " + os.path.abspath(outputFile) + " " + os.path.abspath(outputFile)
  runCmd(cmd, logger)

  logger.info("done.")

def report(logger, configFile, outputFilePrefix, block, dbVersion):
  check_file_exists(configFile)

  rmdScript = os.path.join( os.path.dirname(__file__), "report.rmd")
  if not os.path.exists(rmdScript):
    raise Exception("Cannot find rmdScript %s" % rmdScript)

  if os.path.isfile(dbVersion):
    chromInfo_file = dbVersion
  elif dbVersion == "hg19":
    chromInfo_file = os.path.join(os.path.dirname(__file__), "data", "chromInfo_hg19.txt")
  elif dbVersion == 'hg38':
    chromInfo_file = os.path.join(os.path.dirname(__file__), "data", "chromInfo_hg38.txt")
  else:
    raise Exception("I don't understand dbVersion: %s" % dbVersion)

  level_chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
  level_mut = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']

  logger.info("Reading chromosome length from: %s ..." % chromInfo_file)
  catItems = []# List<CategoryItem>
  chromLengthMap = {}
  with open(chromInfo_file, "rt") as fin:
    for line in fin:
      if line.startswith("#"):
        continue

      parts = line.rstrip().split('\t')
      chrom = parts[0]
      if not chrom in level_chr:
        continue

      chromLength = int(parts[1])
      chromLengthMap[chrom] = chromLength

      chromStart = 0
      while chromStart < chromLength:
        chromEnd = min(chromStart + block, chromLength)
        catItems.append(CategoryItem(chrom, chromStart, chromEnd, None))
        chromStart = chromStart + block

  logger.info("Preparing dinucleotide files ...")
  targetFolder = os.path.dirname(outputFilePrefix)
  targetConfigFile = os.path.join(targetFolder, "dinucleotide_file.list")
  with open(targetConfigFile, "wt") as fout:
    with open(configFile, "rt") as fin:
      for line in fin:
        parts = line.rstrip().split('\t')
        dinuFile = parts[0]
        dinuName = parts[1]
        dinuGroup = parts[2]
        targetDinuFile = os.path.join(targetFolder, "%s_block%d.txt" % (os.path.basename(dinuFile), block))
        fout.write("%s\t%s\t%s\n" % (targetDinuFile, dinuName, dinuGroup))

        if os.path.isfile(targetDinuFile):
          continue

        idxFile = dinuFile + ".tbi"

        check_file_exists(dinuFile)
        check_file_exists(idxFile)

        logger.info("Preparing %s ..." % dinuFile)
        with open(targetDinuFile, "wt") as fdinu:
          tb = tabix.open(dinuFile)

          count = 0
          lastChrom = ""
          for catItem in catItems:
            count = count + 1
            if count % 10000 == 0:
              logger.info("%d / %d" % (count, len(catItems)))

            if catItem.reference_name != lastChrom:
              if lastChrom in chromLengthMap:
                fdinu.write("%s\t%d\t%d\t%d\n" % (lastChrom, chromLengthMap[lastChrom] - 1, chromLengthMap[lastChrom], 0))
              lastChrom = catItem.reference_name
            
            #logger.info("Processing %s:%d-%d ..." %(catItem.reference_name, catItem.reference_start, catItem.reference_end))
            tbiter = tb.query(catItem.reference_name, catItem.reference_start, catItem.reference_end)
            records = [record for record in tbiter]
            totalCount = 0
            for record in records:
              dinucleotide = record[3]
              diCount = int(record[4])
              if diCount == 0:
                diCount = 1
              totalCount += diCount
            
            if totalCount > 0:
              fdinu.write("%s\t%d\t%d\t%d\n" % (catItem.reference_name, catItem.reference_start, catItem.reference_end, totalCount))
          fdinu.write("%s\t%d\t%d\t%d\n" % (lastChrom, chromLengthMap[lastChrom] - 1, chromLengthMap[lastChrom], 0))

  targetRmd = outputFilePrefix + ".rmd"
  shutil.copyfile(rmdScript, targetRmd)

  cmd = "R -e \"setwd('%s');library(knitr);rmarkdown::render('%s',params=list(db='%s'));\"" % (os.path.dirname(targetRmd), os.path.basename(targetRmd), dbVersion)
  runCmd(cmd, logger)
