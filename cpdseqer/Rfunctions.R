# Added with multiQC functions relative to GitHub version 05/20/20
# singleSample(): When only one (group of) sample is specified, performs Chisq test for overall 5 categories and one-vs-other test for four DINUCs.

# INPUT smp: either file name or read-in data.frame for single sample. (When reg is specified, smp must be preprocessed to contain counts for the corresponding region.)
# INPUT gn: string for genome: hg38,hg19,saccer3
# INPUT reg: specifying the single region.
### choice 1 can be a keyword {'intron','exon','promoter'} etc.
### choice 2 can be a underscore-separated vector of four big integers (or proportions), freshly counted from corresponding CUSTOM region from refGenome.
# INPUT cnt5file='gn_reg_cnt5.tsv': one data file pre-storing count vectors for common regions of a few genomes we accommodate.
# OUTPUT res$res: 5x6 test result table in character matrix.
# OUTPUT res$unit: inform on the unit of sample counts used in statistical tests & presentation.
singleSample <- function(smp,gn,reg=NULL,cntType='rCnt',cnt5file='gn_reg_cnt5.tsv0',nDigit=4,unit2=c(1e6,1e3)) { #nDigit: number of significant digits.
	if (!require(reshape2)) stop('Please install library reshape2 for dcast function')
	gn <- tolower(gn)
	DINUC4 <- c('TT','TC','CC','CT')
	if (is.null(reg)) { # Simplest case: no specification of region
		p4gn <- switch(gn,
			hg38=c(0.0986,0.0595,0.0515,0.0700), #0.7204
			hg19=c(0.0980,0.0594,0.0521,0.0700), #0.7206
			saccer3=c(0.108,0.062,0.039,0.058), #0.733
		) # Limited to considered few refGenomes.
	} else { # parse out p4 vector when reg is supplied.
		p4gn <- parseReg(reg,gn,cnt5file=cnt5file)[1:4]
	} # expected probablility vector for five categories.
  if (!is.null(p4gn)) { # Somehow we might fail retrieving expected p-vector
    p5gn <- c(p4gn,1-sum(p4gn))
    names(p5gn) <- c(DINUC4,'others')
		if (!grepl(',',smp)) {
			smpCnt <- collapse.cntFile(smp,cntType,chrs=NULL)
		} else {
			smpCnt <- collapseGrp2Smp(smp)
		}
		if (!is.null(unit2)) {
			res <- convert.cntUnit(smpCnt,unit2)
			smpCnt <- res$cnt
			unit <- res$unit
		} else {
			unit <- 1
		}
    ########### Overall test spanning 5 categories ###########
		all.chisq.res <- chisq.test(smpCnt,p=p5gn)
    ########### Four iterations of ith vs. others ############
		oneVSother <- matrix(NA,nr=4,nc=6,dimnames=list(DINUC4,c('exp','obs','p.exact','p.exact.adj','p.chisq','p.chisq.adj')))
		for (di in DINUC4) {
			x <- smpCnt[di]
			n <- x+smpCnt['others']
			oneVSother[di,'exp'] <- exp <- signif(p5gn[di]/(p5gn[di]+p5gn['others']),nDigit)
			oneVSother[di,'obs'] <- obs <- signif(x/n,nDigit)
			p.binom <- pbinom(x-1,n,exp,lower.tail=F)
			oneVSother[di,'p.exact'] <- signif(max(p.binom,1e-15),nDigit)
			p.chisq <- chisq.test(c(x,smpCnt['others']),p=c(exp,1-exp))$p.val
			oneVSother[di,'p.chisq'] <- signif(max(p.chisq,1e-15),nDigit)
			oneVSother[di,'p.exact.adj'] <- min(signif(p.binom*length(DINUC4),nDigit),1)
			oneVSother[di,'p.chisq.adj'] <- min(signif(p.chisq*length(DINUC4),nDigit),1)
		}
		expCnt <- floor(all.chisq.res$exp)
		smpCnt <- floor(smpCnt)
		expCnt <- paste(expCnt,collapse=':')
		smpCnt <- paste(smpCnt,collapse=':')
		all.row <- c(expCnt,smpCnt,'---','---',signif(all.chisq.res$p.val,nDigit),'---')
		res <- rbind(oneVSother,Overall=all.row)
  } else {
		warning('Cannot retrieve reference dinuc distribution vector for ',gn,' ',reg,'!!!')
		res <- matrix(NA,nr=5,nc=5,dimnames=list(DINUC4,c('exp','obs','95%CI','p.binom','p.chisq')))
  }
	return(list(res=res,unit=unit))
}

# twoReg(): Compare dinuc counts between two region sets in a single sample while taking into account reference dinuc distribution of the two regions in refGenome.
# INPUT smp1: filename for region1-limited sample cnt file.
# INPUT smp2: filename for region2-limited sample cnt file.
# INPUT reg1 & reg2: to inform reference p5 vector. A region keyword or a count5 string (comma-separated).
# NOTE smp1 & smp2: when supplied as multi-samples, smp1 & smp2 should be two groups of equal lengths. Haven't enforced such a pre-check, but it should be kept in mind.
twoReg <- function(smp1,smp2,reg1,reg2,gn,cntType='rCnt',cnt5file='gn_reg_cnt5.tsv',nDigit=4,unit2=c(1e6,1e3),DINUC4=c('TT','TC','CC','CT')) {
	if (!grepl(',',smp1)) {
		smpCnt1 <- collapse.cntFile(smp1,cntType,chrs=NULL)
	} else {
		smpCnt1 <- collapseGrp2Smp(smp1)
	}
	if (!grepl(',',smp2)) {
		smpCnt2 <- collapse.cntFile(smp2,cntType,chrs=NULL)
	} else {
		smpCnt2 <- collapseGrp2Smp(smp2)
	}
	if (!is.null(unit2)) {
		res <- convert.cntUnit(rbind(smpCnt1,smpCnt2),unit2)
		smpCnt1 <- res$cnt[1,]
		smpCnt2 <- res$cnt[2,]
		unit <- res$unit
	} else {
		unit <- 1
	}
	ref5.1 <- parseReg(reg1,gn,scaleNow=F,cnt5file=cnt5file)
	ref5.2 <- parseReg(reg2,gn,scaleNow=F,cnt5file=cnt5file)
	res <- matrix(NA,nr=5,nc=7,dimnames=list(c(DINUC4,'Overall'),c('exp','obs','p.exact','p.exact.adj','p.chisq','p.chisq.adj','1stVS2nd')))
	chisq.stats <- numeric(length(DINUC4))
	names(chisq.stats) <- DINUC4
	for (di in DINUC4) {
		nRef.di <- c(ref5.1[di],ref5.2[di])
		p.di <- nRef.di/sum(nRef.di)
		obs.di <- c(smpCnt1[di],smpCnt2[di])
		chisq.res <- chisq.test(x=obs.di,p=p.di)
		p.chisq <- chisq.res$p.val
		chisq.stats[di] <- chisq.res$stat
		res[di,'exp'] <- paste(floor(nRef.di/unit),collapse=':') # Showing ratio (a:b)
		res[di,'obs'] <- paste(floor(obs.di),collapse=':') # Showing ratio (a:b)
		p.binom <- pbinom(obs.di[1]-1,sum(obs.di),p.di[1],lower.tail=F)
		res[di,'p.chisq'] <- p.chisq <- signif(max(p.chisq,1e-15),nDigit)
		res[di,'p.exact'] <- p.binom <- signif(max(p.binom,1e-15),nDigit)
		res[di,'p.exact.adj'] <- min(signif(p.binom*length(DINUC4),nDigit),1)
		res[di,'p.chisq.adj'] <- min(signif(p.chisq*length(DINUC4),nDigit),1)
		res[di,'1stVS2nd'] <- ifelse((obs.di[1]-1)/sum(obs.di)>p.di[1],'>','<')
	}
	p.exact.all <- sumPbyFisher(as.numeric(res[DINUC4,'p.exact']))
	res['Overall','p.exact'] <- signif(max(p.exact.all,1e-15),nDigit)
	p.chi2.all <- sumPbyChi2(chisq.stats)
	res['Overall','p.chisq'] <- signif(max(p.chi2.all,1e-15),nDigit)
	res['Overall',setdiff(colnames(res),c('p.exact','p.chisq'))] <- '---'
	return(list(res=res,unit=unit))
}
# twoGrp(): Given two sample file names (.cnt), return 5*7 statistical test result table.
# INPUT smp1: either file name or read-in data.frame for sample 1.
# INPUT smp2: either file name or read-in data.frame for sample 2.
# OUTPUT res$res: 5x7 test result table.
### Five rows (DINUC4+Overall), 7 columns (exp,obs,p.exact,p.exact.adj,p.chisq,p.chisq.adj,1stVS2nd)
# NOTE only chisquared test is employed for both overall and one-vs-others.
# NOTE number of digits printed out in kable is adjustable via nDigit
twoGrp <- function(smp1,smp2,cntType='rCnt',nDigit=4,DINUC4=c('TT','TC','CC','CT'),unit2=c(1e6,1e3)) {
	if (!grepl(',',smp1)) {
		smpCnt1 <- collapse.cntFile(smp1,cntType,chrs=NULL)
	} else {
		smpCnt1 <- collapseGrp2Smp(smp1)
	}
	if (!grepl(',',smp2)) {
		smpCnt2 <- collapse.cntFile(smp2,cntType,chrs=NULL)
	} else {
		smpCnt2 <- collapseGrp2Smp(smp2)
	}
	if (!is.null(unit2)) {
		res <- convert.cntUnit(rbind(smpCnt1,smpCnt2),unit2)
		smpCnt1 <- res$cnt[1,]
		smpCnt2 <- res$cnt[2,]
		unit <- res$unit
	} else {
		unit <- 1
	}
	res <- matrix(NA,nr=5,nc=7,dimnames=list(c(DINUC4,'Overall'),c('exp','obs','p.exact','p.exact.adj','p.chisq','p.chisq.adj','1stVS2nd')))
	for (di in DINUC4) {
		N2.1 <- c(smpCnt1[di],smpCnt1['others']) # two numbers of grp1
		N2.2 <- c(smpCnt2[di],smpCnt2['others'])# two numbers of grp2
		chisq.res <- chisq.test(rbind(N2.1,N2.2))
		exp1 <- floor(chisq.res$expected[1,1])
		exp2 <- floor(chisq.res$expected[2,1])
		res[di,'exp'] <- paste(exp1,exp2,sep=', ')
		res[di,'obs'] <- paste(floor(N2.1[1]),floor(N2.2[1]),sep=', ') # Observed di cnts in millions
		p.chisq <- max(chisq.res$p.val,1e-15) # When too little p (say 0) occurs, upgrade to 1e-15.
		res[di,'p.chisq'] <- signif(p.chisq,nDigit)
		res[di,'p.chisq.adj'] <- min(signif(p.chisq*length(DINUC4),nDigit),1)
		m=smpCnt1[di]+smpCnt2[di]
		n=smpCnt1['others']+smpCnt2['others']
		k=sum(N2.1)
		x=smpCnt1[di]
		p.hyper1 <- phyper(x-1,m,n,k,lower.tail=F)
		p.hyper2 <- phyper(x,m,n,k)
		p.hyper <- max(min(p.hyper1,p.hyper2),1e-15)
		res[di,'p.exact'] <- signif(p.hyper,nDigit)
		res[di,'p.exact.adj'] <- min(signif(p.hyper*length(DINUC4),nDigit),1)
		res[di,'1stVS2nd'] <- ifelse(N2.1[1]>=chisq.res$expected[1,1],'>','<')
	}
	all.chisq.res <- chisq.test(rbind(smpCnt1,smpCnt2))
	obs1 <- paste(floor(smpCnt1),collapse=':')
	obs2 <- paste(floor(smpCnt2),collapse=':')
	exp.M <- floor(all.chisq.res$exp)
	exp1 <- paste(exp.M[1,],collapse=':')
	exp2 <- paste(exp.M[2,],collapse=':')
	res['Overall','exp'] <- paste(exp1,exp2,collapse='; ')
	res['Overall','obs'] <- paste(obs1,obs2,collapse='; ')
	res['Overall','p.exact'] <- res['Overall','p.exact.adj'] <- res['Overall','p.chisq.adj'] <- '---'
	res['Overall','p.chisq'] <- signif(max(all.chisq.res$p.val,1e-15),nDigit)
	res['Overall','1stVS2nd'] <- ifelse(smpCnt1[1]>=all.chisq.res$exp[1,1],'>','<')
	return(list(res=res,unit=unit))
}

# collapse.cntFile(): Given one cnt filename, derive smpCnt (five cnts) vector.
# (interim) UPDATE 5/10/20: Enforce 3rd column named "ReadCount"
# UPDATE 5/9/20: add option cntType with default value of rCount.
# INPUT smp: sample file name (.cnt) or a data frame
collapse.cntFile <- function(smp,cntType=c('rCnt','sCnt')[1],chrs=NULL,DINUC4=c('TT','TC','CC','CT')) {
	if (!require(reshape2))
		stop('Please install library reshape2 for dcast function')
	if (length(as.vector(smp))==1)
		smp <- read.delim(smp,as.is=T)
	colnames(smp)[3] <- 'ReadCount' # UPDATE 5/10/20 to accomodate old-versioned CNT file.
	if (!is.null(chrs))
		smp <- smp[smp[,'Chromosome']%in%chrs,]
	dinuc16 <- genPerm2()
	smp <- smp[smp[,'Dinucleotide']%in%dinuc16,]
	if (cntType=='rCnt') {
		smp <- reshape2::dcast(smp,Dinucleotide~Chromosome,value.var='ReadCount')
	} else if (cntType=='sCnt') {
		smp <- reshape2::dcast(smp,Dinucleotide~Chromosome,value.var='SiteCount')
	} else {
		stop('Dinucleotide count type must be either rCnt or sCnt!')
	}
	smpMat <- as.matrix(smp[,-1])
	rownames(smpMat) <- smp$Dinucleotide
	colnames(smpMat) <- colnames(smp)[-1]
	smpCnt0 <- apply(smpMat,1,function(x) sum(as.numeric(x),na.rm=T))
	smpCnt.4 <- smpCnt0[DINUC4]
	smpCnt.5th <- sum(as.numeric(smpCnt0[setdiff(names(smpCnt0),DINUC4)]))
	smpCnt <- c(smpCnt.4,others=smpCnt.5th)
}
# collapse16.cntFile(): Given one cnt filename, sum up counts for all 16 dinucleotide types.
# INPUT smp: sample file name (.cnt) or a data frame
collapse16.cntFile <- function(smp,cntType=c('rCnt','sCnt')[1],chrs=NULL,DINUC4=c('TT','TC','CC','CT')) {
  if (!require(reshape2))
    stop('Please install library reshape2 for dcast function')
  if (length(as.vector(smp))==1)
    smp <- read.delim(smp,as.is=T)
	colnames(smp)[3] <- 'ReadCount' # UPDATE 5/10/20 to accomodate old-versioned CNT file.
  if (!is.null(chrs))
    smp <- smp[smp[,'Chromosome']%in%chrs,]
  dinuc16 <- genPerm2()
  smp <- smp[smp[,'Dinucleotide']%in%dinuc16,]
  if (cntType=='rCnt') {
    smp <- reshape2::dcast(smp,Dinucleotide~Chromosome,value.var='ReadCount')
  } else if (cntType=='sCnt') {
    smp <- reshape2::dcast(smp,Dinucleotide~Chromosome,value.var='SiteCount')
  } else {
    stop('Dinucleotide count type must be either rCnt or sCnt!')
  }
  smpMat <- as.matrix(smp[,-1])
  rownames(smpMat) <- smp$Dinucleotide
  colnames(smpMat) <- colnames(smp)[-1]
  smpCnt <- apply(smpMat,1,function(x) sum(as.numeric(x),na.rm=T))
  smpCnt
}

# genPerm2(): # generate permutations of two-tuples.
genPerm2 <- function(vocab=c('A','T','G','C')) {
	comb.2col <- t(combn(vocab,2))
	comb.col2 <- comb.2col[,2:1]# reverse 1st col & 2nd col of comb.2col
	sameRepeat <- cbind(vocab,vocab)
	prototype <- data.frame(rbind(comb.2col,comb.col2,sameRepeat))
	perm2 <- do.call(paste0,prototype)
	perm2
}
# parseReg(): parse reg to retrieve p5 vector. Typically extraction of the first 4 values will immediately follow.
# INPUT reg: either one keyword or five underscore separated frequency values.
# OUTPUT p5gn: a 4-component frequency vector. Program can abort at unexpected input.
# INPUT scaleNow: logic variable instructing if an across-5-dinuc scaling should be performed at the current step. set to TRUE by default.
parseReg <- function(reg,gn='hg38',scaleNow=T,cnt5file='gn_reg_cnt5.tsv',DINUC4=c('TT','TC','CC','CT')) {
	reg <- separated <- unlist(strsplit(reg,','))
	if (length(separated)==5) { # Underscore separated into 5 numeric strings.
		separated <- as.numeric(separated)
		p5gn <- separated
		if (max(p5gn)>1 & scaleNow) p5gn <- separated/sum(as.numeric(separated))
	} else if (length(separated)==1) {
		p5mat <- read.delim(cnt5file,row.names=1)
		if (is.element(paste(gn,reg,sep='_'),rownames(p5mat)))
			p5gn <- as.numeric(p5mat[paste(gn,reg,sep='_'),])
			if (scaleNow) p5gn <- p5gn/sum(p5gn)
	} else {
		stop('Please supply one keyword or five underscore-separated frequency values for one region!!!')
	}
	names(p5gn) <- c(DINUC4,'others')
	p5gn
}

# INPUT chi2.dinucs: chisq test statistics across two regions for each dinuc.
sumPbyChi2 <- function(chi2.dinucs) {
	sumchi2 <- sum(chi2.dinucs)
	dof <- length(chi2.dinucs)
	p <- pchisq(sumchi2,dof,lower.tail=F)
	p
}
# INPUT pvals: individual test p-values across two regions for each dinuc.
sumPbyFisher <- function(pvals,min.th=1e-147) {
	pvals[pvals<min.th] <- min.th
	aggChi2 <- (-2*sum(log(pvals)))
	dof <- 2*length(pvals)
	p <- pchisq(aggChi2,dof,lower.tail=F)
	p
}

# collapseGrp2Smp(): when sample CNT file name implies multiple files are supplied, flattenGrp2Smp() is invoked to collapse a group to a sample.
# NOTE: collapseGrp2Smp() calls on collapse.cntFile() to yield sample-wise 5-tuple count vectors, and then sum up counts by dinuc type.
# INPUT grp: a comma separated string, implying multiple CNT files.
# OUTPUT CNT5: a 5-tuple of counts, named DINUC4 + others.
collapseGrp2Smp <- function(grp) {
	samples <- unlist(strsplit(grp,','))
	cnt.samples <- matrix(nr=length(samples),nc=5)
	for (i in 1:length(samples)) {
		sample <- samples[i]
		cnt5 <- collapse.cntFile(sample)
		cnt.samples[i,] <- cnt5
	}
	CNT5 <- apply(cnt.samples,2,function(x) sum(x,na.rm=T))
	names(CNT5) <- names(cnt5)
	CNT5
}
# convert.cntUnit(): Convert dinuc count to bigger unit, say a million.
# INPUT cnt0: raw counts, most likely as a 5-value vector.
# INPUT unit2: two candidate count units, one bigger (say 1m) one smaller (say 1k). 
# OUTPUT res: list of two components
###### $cnt - counts converted to bigger unit.
###### $unit - the actual unit used in converted counts.
convert.cntUnit <- function(cnt0,unit2=c(1e6,1e3)) {
	uH <- max(unit2)
	uL <- min(unit2)
	cnt <- floor(cnt0/uH)
	unit <- uH
	if (min(cnt)<1) {
		cnt <- floor(cnt0/uL)
		unit <- uL
	}
	res <- list(cnt=cnt,unit=unit)
	res
}


################################################################################################
############################ QC functions below ################################################
################################################################################################


# efficiency(): Calculate efficiency of a single sample (CNT file).
# smp: sample CNT file or a 6-columned bed data frame.
# cntType: allow toggling between rCnt (default) & sCnt.
efficiency <- function(smp,cntType=c('rCnt','sCnt')[1],nDigit=3) {
	DINUC4=c('TT','TC','CC','CT')
	if (length(as.vector(smp))==1) {
		cnt5 <- collapse.cntFile(smp,cntType)
	} else { # smp is a bed data.frame
		cnt16 <- bed2cnt(smp,cntType)
		cnt5 <- c(cnt16[DINUC4],sum(cnt16[setdiff(names(cnt16),DINUC4)],na.rm=T))		
	}
	dinuc4 <- sum(cnt5[DINUC4],na.rm=T)
	eff <- sum(dinuc4,na.rm=T)/sum(cnt5,na.rm=T)
	eff <- signif(eff,nDigit)
	eff 
}

# contrast(): Calculate contrast for a single sample (CNT file).
# smp: sample CNT file or a 6-columned bed data frame.
# cntType: allow toggling between rCnt (default) & sCnt.
contrast <- function(smp,cntType=c('rCnt','sCnt')[1],nDigit=3) {
	if (length(as.vector(smp))==1) {
		cnt16 <- collapse16.cntFile(smp,cntType)
	} else {
		cnt16 <- bed2cnt(smp,cntType)
	}
	cntTT <- as.numeric(cnt16['TT']) # Removed name of value
	cnt0 <- as.numeric(cnt16['AA'])
	if (cnt0>0) {
		cont <- cntTT/cnt0 
	} else {
		cont <- NA
	}
	cont <- signif(cont,nDigit)
	cont
}

# symmetry(): to derive 7*2 values for subsequent symmetry plotting.
# INPUT bgz: BGZ file name.
# OUTPUT sym7: counts of TT,TC,CC,CT,DINUC4, efficiency, and contrast. One series for + strand, the other for - strand. 
symmetry <- function(bgz,cntType=c('rCnt','sCnt')[1],nDigit=3) {
	if (!require(data.table)) {
		warning('No data.table available; time-consuming to read in BGZ files!')
		bed <- read.delim(bgz,header=F,as.is=T)
	}	else {
		bed <- fread(cmd=paste0(c('zcat',bgz),collapse=' '),header=F,data.table=F)
	}
	colnames(bed) <- c('Chromosome','start0','End','Dinucleotide','Count','strand') #start0 - open start point.
	bedP <- bed[bed[,'strand']=='+',]
	bedN <- bed[bed[,'strand']=='-',]
	cnt16P <- bed2cnt(bedP,cntType)
	cnt16N <- bed2cnt(bedN,cntType)
	DINUC4=c('TT','TC','CC','CT')
	cnt4 <- rbind(cnt16P[DINUC4],cnt16N[DINUC4])
 	eff <- sapply(list(bedP,bedN),efficiency,cntType,nDigit)
	cont <- sapply(list(bedP,bedN),contrast,cntType,nDigit)
	sym7 <- cbind(cnt4,rowSums(cnt4),eff,cont)
	colnames(sym7) <- c(DINUC4,'DINUC4','Efficiency','Contrast')
	rownames(sym7) <- c('Forward','Reverse')
	list(sym7=sym7,bed=bed)	
}

# bed2cnt(): derive counts for 16 dinucs from a bgz-formatted bed data.frame.
### Accompolishes second task of bam2dinucleotide (.count) and collapse.cntFile() task.
# INPUT bed: data.frame of 6 cols - c('Chromosome','start0','End','Dinucleotide','Count','strand') #start0 - open start point.
# INPUT chrs can be set as paste0('chr',c(1:22,c('X','Y','M')))
bed2cnt <- function(bed,cntType=c('rCnt','sCnt')[1],chrs=NULL) {
  #library(reshape2)
	if (!is.null(chrs))
	  bed <- bed[bed[,1]%in%chrs,]
	dinuc16 <- genPerm2()
	bed <- bed[bed[,'Dinucleotide']%in%dinuc16,]
	if (cntType=='rCnt') {
	  cnt <- tapply(bed$Count,bed[,'Dinucleotide'],sum) #by is also OK
	} else if (cntType=='sCnt') {
		cnt <- tapply(bed$Count,bed[,'Dinucleotide'],length) #by is also OK
	} else {
		stop('Dinucleotide count type must be either rCnt or sCnt!')
	}
	cnt <- cnt[dinuc16]
	cnt
}
#smp <- commandArgs(T)
#symRes <- symmetry('case1.bgz')
#bed <- symRes$bed
#sym7 <- symRes$sym7

# plot_EfforCont(): Plot left-and-right two panels for Eff & Cont, respectively.
### NOTE: The following four reference intervals are determined from ten human CPD-seq samples as of May 2020.
EffR.r=c(0.576,0.813) #EffRange.rCnt
EffR.s=c(0.575,0.805) #EffRange.sCnt
ContR.r=c(7.73,23.8) #ContRange.rCnt
ContR.s=c(7.47,23.2) #ContRange.sCnt
plot_EffCont <- function(eff2,cont2,sample='Sample',EffRange.rCnt=EffR.r,EffRange.sCnt=EffR.s,ContRange.rCnt=ContR.r,ContRange.sCnt=ContR.s) {
	names(eff2) <- names(cont2) <- c('by_ReadCount','by_SiteCount')
	layout(matrix(1:2,nr=1))
	plot_DotByStem(eff2,sample=sample,index='Efficiency',range.rCnt=EffRange.rCnt,range.sCnt=EffRange.sCnt)
	plot_DotByStem(cont2,sample=sample,index='Contrast',range.rCnt=ContRange.rCnt,range.sCnt=ContRange.sCnt)
}

# plot_DotByStem(): plot stem-and-dot figure for one index, either Efficiency or Contrast.
# INPUT indexVal: two values of the index, named c('by_ReadCount','by_SiteCount')
# INPUT range.rCnt (byReadCount) of Efficiency out of 10 human samples: c(0.576,0.813)
# INPUT range.rCnt (byReadCount) of Contrast out of 10 human samples: c(7.7,23.8) 

plot_DotByStem <- function(indexVal,range.rCnt=c(0.576,0.813),range.sCnt=c(0.576,0.813), sample='Sample',index='Efficiency') {
	par(mar=c(6,8,4,2),cex.main=1.5,cex.axis=1.5,cex.lab=1.5,bty='l')#mar=c(6,8,4,2),
	xInterval <- switch(index,
		Efficiency=c(0,1),
		Contrast=c(0,30)
	)
	x.offset <- switch(index,
		Efficiency=-0.25,
		Contrast=-5
	)
		plot(xInterval,c(0,3),type='n',xlab=index,ylab='',main=sample,axes=F)
	axis(side=1)
	#axis(side=2,at=c(0:3),labels=c('',rev(names(indexVal)),''),las=1)
	text(cex=1.5, x=x.offset, y=c(1,2)-0.25, rev(names(indexVal)), xpd=TRUE, srt=45)
	lines(range.rCnt,c(2,2),lwd=4) # Set by_ReadCount higher, at y=2
	points(indexVal[1],2,pch=19,col='red',cex=2)
	lines(range.sCnt,c(1,1),lwd=2) # Set by_SiteCount lower, at y=1
	points(indexVal[2],1,pch=4,col='blue',cex=2,lwd=4)
}

# plot_sym7(): Given 2-by-7 indices of symmetry, plot seven horizontal bars to signify symmetry situation.
plot_sym7 <- function(sym7,sample='') {
	par(mar=c(6,8,4,10),fin=c(12,6),cex.axis=1.5,cex.main=1.5) # UPDATE 05/21/20 fin
	sym7.normed <- t(t(sym7)/colSums(sym7))
	barplot(sym7.normed[,ncol(sym7.normed):1],horiz=T,las=1,main=paste(sample,'Symmetry'),legend.text=T,
		args.legend=list(x='right',pt.cex=2,cex=1.5,bty='n',inset=c(-0.8,0)))
	abline(v=0.5,lty='dashed',col='red',lwd=2)
}

# plot_rcDistrib(): Take the bed data frame as primary input, plot a groupped barplot for three levels of readCounts.
# INPUT bed0: secondary output of symmetry().
# INPUT pts: default to c(0,5,10) points where frequency of higher readCounts is visualized. 
# INPUT sample: sample name, to be prefixed to title (main) of figure.
# OUTPUT Freq: Five-by-three count frequncy table. Top four rows for TT, TC, CC, and CT, bottom row for All. 
######## 5 Rows: TT, TC, CC, CT, and All
######## 3 Columns: >0, >5, >10
plot_rcDistrib <- function(bed0,sample='',pts=c(0,5,10)) {
	rc0 <- bed0[,5]
  DINUC4=c('TT','TC','CC','CT')
	bed <- bed0[bed0[,4]%in%DINUC4,]
  dinucs <- bed[,4]
  rc <- bed[,5]
  RC4 <- split(rc,factor(dinucs))[DINUC4]
	#Freq <- matrix(nr=length(DINUC4)+1,nc=length(pts),dimnames=list(c(DINUC4,'All'),paste('Reads',pts,sep='>')))
	Freq <- sapply(pts,function(x) sapply(RC4,exceedingCnt,x))
	Freq.all <- exceedingCnt(rc0,pts)	
	colnames(Freq) <- paste('Reads',pts,sep='>')
	par(mar=c(6,8,6,2),fin=c(12,6),cex.main=1.5,cex.axis=1.5) # UPDATE 05/21/20 fin
	barplot(Freq,beside=T,log='y',las=1,main=paste(sample,'\nRead Count Frequency'),legend.text=T,
		args.legend=list(x='topright',pt.cex=1.5,cex=1.2))
	mtext(side=2,line=5,text='Frequency',cex=2)
  #tiff(paste0('rcDistrib_truncate',truncate,'.tif'),width=2048,height=1600)
  #dev.off()
  rbind(Freq,All=Freq.all) 
}

# exceedingCnt(): How many values of rc exceed (>) the specified value of pt? 
# UPDATE 5/15/20: Allow pt to be flexible, scalar or vector
exceedingCnt <- function(rc,pt) {
	freq <- sapply(pt, function(x)  sum(rc>x,na.rm=T))
	freq
}
# plotM_rcDistrib(): Given a list of bgz bed data frames, plot barplots of readCounts for multiple samples, at given cut-off points.
# INPUT bedS: a list of K components, corresponding to K samples in a batch.
# dinuc takes value from {'TT','TC','CC','CT'} or any combination
plotM_rcDistrib <- function(bedS,exp='Experiment',dinuc=c('TT','TC','CC','CT'),pts=c(0,5,10)) {
	K <- length(bedS)
	bedS <- lapply(bedS,function(x,scope) x[x[,4]%in%scope,],dinuc)
	rcS <- lapply(bedS,function(x) x[,5])	
	Freq <- sapply(rcS,exceedingCnt,pts) # count matrix: pts by sample
	Freq <- t(Freq) # transpose to sample by pts
	colnames(Freq) <- paste('Reads',pts,sep='>')
	Freq <- Freq[order(-Freq[,1]),]
  par(mar=c(6,8,6,2),fin=c(12,6),cex.main=1.5,cex.axis=1.5) #mar=c(6,8,6,2),
  barplot(Freq,beside=T,log='y',las=1,main=paste(exp,'\nRead Count Frequency'),legend.text=T,
    args.legend=list(x='topright',bty='n',pt.cex=1.2))
  mtext(side=2,line=5,text='Frequency',cex=2)
	Freq
}
# plotM_EffCont(): Calculate Efficiency & Contrast for multiple samples and plotted them in left-and-right panels.
# INPUT effs: Efficiency values for K samples.
# INPUT conts: Contrast values for K samples.
plotM_EffCont <- function(effs,conts,range.eff,range.cont,exp='Experiment') {
	layout(matrix(1:2,nr=1))
	plot_particles(effs,'Efficiency',range.eff,exp)
	plot_particles(conts,'Contrast',range.cont,exp)	
}

# plot_particles(): Spread out K Eff/Cont values vertically as particles. Customized boxplot.
# INPUT indexVal: K index values, must be named with sample names.
plot_particles <- function(indexVal,index='Efficiency',indexRange=c(0.576,0.813),Experiment='') {
	indexVal <- sort(indexVal[!is.na(indexVal)])
	K <- length(indexVal)
	par(mar=c(4,9,4,2),cex.main=1.5,cex.axis=1.2,cex.lab=1.5)
	plot(range(indexVal),c(-0.5,K),type='n',axes=F,ylab='',xlab=index,main=Experiment)
	num5 <- fivenum(indexVal)
	q25 <- num5[2]
	q75 <- num5[4]
	iqrVal <- IQR(indexVal)
	lb <- q25-1.5*iqrVal
	ub <- q75+1.5*iqrVal
	outliers <- indexVal<lb | indexVal>ub
	for (i in 1:K) {
		if (!outliers[i]) {
			points(indexVal[i],i,pch=19,cex=1.5)
		} else {
			points(indexVal[i],i,pch='x',cex=2,col='red')
		}
	}	
	lines(indexRange,c(0,0),col='gray',lwd=2)
	text(median(indexVal),-0.4,'reference',cex=1.2,col='gray')
	axis(side=1)#,at=seq(from=min(indexVal),to=max(indexVal),len=5))
	axis(side=2,at=1:K,labels=names(indexVal),las=1)
}
