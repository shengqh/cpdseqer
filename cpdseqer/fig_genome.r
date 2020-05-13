inputFile='genome.config.txt'
chromInfoFile='chromInfo.txt'
block='100000'
outfilePrefix='genome'
setwd('/gpfs23/scratch/cqs/shengq2/guoyan/20191104_cpd_analysis/fig_genome')

plot.pos = TRUE
plot.mut = TRUE
ctrl.default = FALSE

library(data.table) ## fread
library(dplyr)
library(ggplot2)
library(reshape2) ## acast

samples <- read.delim(inputFile, sep = '\t', header = T, stringsAsFactors = F)

if(nrow(samples) == 0) {
  stop(paste0("No data in ", inputFile))
}

if(!all(samples$Group %in% c(0, 1))){
  stop("Group should be either 0 for control or 1 for case in dinucleotide_file.list")
}

controlSamples<-samples$Sample[samples$Group==0]
caseSamples<-samples$Sample[samples$Group==1]
doBinomialTest<-(length(controlSamples) > 0) & (length(controlSamples > 0))

level.chr <- c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY')
level.mut <- c('AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT')

chrom_levels<-rev(level.chr)

info <- read.delim("chromInfo.txt", sep = "\t", header = TRUE, stringsAsFactors = F)

colnames(info)<-c(paste0("V", c(1:ncol(info))))
info <- info[info$V1 %in% level.chr, c("V1", "V2")]
info <- data.frame(chrom = factor(info$V1, levels = chrom_levels),
                   chromStart = rep(0, nrow(info)),
                   chromEnd = info$V2)

ctype = "Site"
for (ctype in c("Site", "Read")){
  columnName = paste0(ctype, "Count")
  curPrefix = paste0(outfilePrefix, "_", columnName)
  pdf(paste0(curPrefix, ".pdf"), width=7, height=7, onefile = T)
  ## cases plot
  i<-1
  cnt.all <- NULL
  for (i in 1:nrow(samples)) {
    sampleFile <- samples$DinuFile[i]
    countFile <- samples$CountFile[i]
    sampleName <- samples$Sample[i]
    sampleGroup <- samples$Group[i]

    case <- fread(sampleFile, sep = "\t", header = T, data.table = F, stringsAsFactors = F)
    case <- case[case$Chrom %in% level.chr,]
    case$Chrom<-factor(case$Chrom, levels=chrom_levels)

    seg <- ggplot(data=info,aes(y = chrom, x = chromStart))+
      geom_segment(aes(y = chrom, yend = chrom, x = chromStart, xend = chromEnd),
                  lineend = "round", color = "Gainsboro", size = 5)+
      geom_segment(data = case,
                  aes_string(y = "Chrom", yend = "Chrom", x = "Start", xend = "End", color=columnName),
                  lineend = "butt", size = 5)+
      scale_y_discrete(drop=FALSE)+
      scale_color_gradient(high="DimGray",low="Gainsboro", name = paste0("Count per ", as.numeric(block) / 1000, " KB")) +
      scale_x_continuous(limits = c(0, 250e6),
                        expand = c(0, 2e6), ## expand dist = (maxlimit-minlimit)*a + b
                        breaks = seq(0, 250e6, by = 50e6),
                        labels = paste0(seq(0, 250, by = 50), "Mb"),
                        position = "top")+
      theme_minimal()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = margin(0.1, 10, 0.1, 0.1, "pt"))+
      labs(x = "Genomic positions", y = NULL, title=paste0("Sample: ", sampleName))
    print(seg)

    df<-fread(countFile, sep="\t", header=T, data.table=F, stringsAsFactors = F)
    df<-df[df$Chromosome %in% level.chr,]
    df<-df[df$Dinucleotide %in% level.mut,]

    df<-df[,c("Chromosome", "Dinucleotide", columnName)]
    df<-dcast(df, Dinucleotide~Chromosome,value.var=columnName)
    rownames(df)<-df$Dinucleotide
    df<-df[,c(2:ncol(df))]
    
    df <- cbind(df, total = apply(df, 1, sum))

    curSample = data.frame("mut" = rownames(df), "count"=df$total)
    colnames(curSample)<-c("mut", sampleName)

    if (is.null(cnt.all)){
      cnt.all = curSample
    }else{
      cnt.all = merge(cnt.all, curSample, id="mut")
    }
    
    totalCount<-sum(df$total)
    df <- apply(df, 2, function(x) x/sum(x))
    
    lab <- rownames(df)
    lab.pos <- cumsum(df[,1]) - df[,1]/2

    df <- melt(df)
    colnames(df) <- c("mut", "chrom", "value")

    bar <- ggplot(df, aes(x = factor(chrom, levels = rev(c(level.chr, "total"))),
                          y = value,
                          fill = factor(mut, levels = rev(level.mut))))+
      geom_bar(color = "black", stat = "identity", width = 0.8)+
      annotate(geom = "text", x = Inf, y = lab.pos, label = lab, size = 4) +
      labs(x = NULL, y = NULL, title = sampleName)+
      scale_y_continuous(expand = c(0.01,0),
                        labels = scales::percent)+
      scale_x_discrete(expand = c(0.035, 0))+
      scale_fill_manual(values = rep(c("white", "Gainsboro"), 8), drop=FALSE)+
      theme_void()+
      theme(axis.text = element_text(size = 10),
            plot.margin = margin(20, 10, 2, 2, unit = "pt"),
            legend.position = "none")+
      coord_flip(clip = "off")
    print(bar)
  }
  dev.off()

  controlSamples<-samples$Sample[samples$Group==0]
  caseSamples<-samples$Sample[samples$Group==1]

  controls<-cnt.all[,colnames(cnt.all) %in% controlSamples,drop=F]
  controlCount<-apply(controls, 1, sum)
  controlCount<-controlCount/sum(controlCount)
  ctrl.cnt<-data.frame('n'=controlCount)
  rownames(ctrl.cnt)<-cnt.all$mut

  cases<-cnt.all[,colnames(cnt.all) %in% caseSamples,drop=F]
  caseCount<-apply(cases, 1, sum)
  case.cnt<-data.frame('n'=caseCount)
  rownames(case.cnt)<-cnt.all$mut

  res.bit <- data.frame(matrix(nrow = nrow(case.cnt), ncol = 0))
  totalCaseCount=sum(case.cnt[, 'n'])
  for (i in 1:nrow(case.cnt)) {
    mut <- rownames(case.cnt)[i]
    bit <- tryCatch(binom.test(case.cnt[mut, 'n'], totalCaseCount, p = ctrl.cnt[mut, 'n']),
                    error = function(e) return(NA))
    res.bit$Mut[i] <- mut
    res.bit$Count[i]<-case.cnt[mut, 'n']
    res.bit$TotalCount[i]<-totalCaseCount
    if(!is.na(bit)[1]){
      res.bit$Proportion[i] <- bit$estimate
      res.bit$Expected[i] <- bit$null.value
      res.bit$pvalue[i] <- bit$p.value
      res.bit$`0.95CI`[i] <- paste0('[', round(bit$conf.int[1], 4), ';', round(bit$conf.int[2], 4), ']')
    }else{
      res.bit$Proportion[i] <- case.cnt[mut, 'n']/totalCaseCount
      res.bit$Expected[i] <- ctrl.cnt[mut, 'n']
      res.bit$pvalue[i] <- NA
      res.bit$`0.95CI`[i] <- NA
    }
  }

  write.table(res.bit, file=paste0(curPrefix, ".txt"), sep="\t", row.names=F)
}