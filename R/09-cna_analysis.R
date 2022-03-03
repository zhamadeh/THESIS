######################################################################
####   Packages   ####
######################################################################

library(rtracklayer)
library(tidyverse)
library(dplyr)
library(plyr)
library(tidyr)
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(GenomicRanges))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(GenomicAlignments))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(DNAcopy))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(doParallel))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(tiydyverse))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(ggplot2))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(cowplot))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(caTools))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(gplots))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(ReorderCluster))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(BiocManager))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(BiocManager::install("AneuFinder"))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(AneuFinder))))
######################################################################
####   Assemble datasets   ####
######################################################################

cutoff=1000
file = "/Users/zeidh/Desktop/binsize_1e+06_PE_CNV.bed.gz"

copyNumberFinder <- function(cutoff, file){
	
	CNV<- import(file)
	cnvPerCell<-data.frame()
	cnvPerCellSummary=data.frame()
	
	for (i in 1:length(CNV)){
		
		## ASSIGN GENOTYPE BASED OFF FILE NAME
		file=strsplit(CNV[i]@listData[[1]]@metadata$trackLine@description,split = " ")[[1]][4]
		
		message("Reading file: ",file," ... ",round((i/length(CNV))*100,2),"%")
		
		## RETRIEVE FIRST LIBRARY AND REMOVE SMALL SEGMENTS < 500Kb and HIGH PLOIDY STATES TO ASSIGN  PLOIDY
		tmp <- as.data.frame(CNV[i])
		tmp$name <- as.factor(tmp$name)
		tmp=tmp[tmp$width>100000,]
		tmp=tmp[tmp$seqnames!="chrY",]
		#HIGH PLOIDY STATES
		tmp <- tmp[tmp$name!="20-somy" & tmp$name!="19-somy" &tmp$name!="18-somy" &tmp$name!="17-somy" & tmp$name!="16-somy" &tmp$name!="15-somy" &tmp$name!="14-somy"&tmp$name!="13-somy" &tmp$name!="12-somy"&tmp$name!="11-somy"&tmp$name!="10-somy"&tmp$name!="zero-inflation",]
		tmp$name <-droplevels(tmp$name)
		
		
		## BUILD PLOIDY TABLE SUMMARY TO CLASSIFY NATIVE STATE AND GAINS & LOSSES
		ploidyTable <- tmp %>% group_by(name) %>% dplyr::summarize(n(),sum=sum(width))
		ploidy = as.numeric(strsplit(as.character(ploidyTable[which.max(ploidyTable$sum),]$name),split = "[-]")[[1]][1])
		if (ploidy==0){
			ploidy=1
		} else if (ploidy >2){
			next
		}
		
		
		
		## REMOVE SMALL SEGMENTS  AFFTER HAVING  ASSIGNED PLOIDY  STATE
		tmp=tmp[tmp$width>cutoff,]
		ploidyTable <- tmp %>% group_by(name) %>% dplyr::summarize(n(),sum=sum(width))
		
		
		## TURN PLOIDY INTO NUMERIC AND CLASSIFY  GAINS/LOSSES
		ploidyTable=  separate(ploidyTable,col = name,sep="-",into = c("name"))
		ploidyTable$name<- as.numeric(ploidyTable$name)
		ploidyTable =  ploidyTable[ploidyTable$name!=ploidy,]
		ploidyTable$type =  ifelse(ploidyTable$name < ploidy, "loss", ifelse((ploidyTable$name > ploidy), "gain", "unclear"))
		## REPEAT FOR WHOLE DF
		tmp=  separate(tmp,col = name,sep="-",into = c("name"))
		tmp$name<- as.numeric(tmp$name)
		tmp =  tmp[tmp$name!=ploidy,]
		tmp$type =  ifelse(tmp$name < ploidy, "loss", ifelse((tmp$name > ploidy), "gain", "unclear"))
		
		
		## QUANTIFY GAINS/LOSSES
		totalGain= sum(filter(ploidyTable,type=="gain")$sum)
		totalGainSeg = sum(filter(ploidyTable,type=="gain")$`n()`)
		totalLoss= sum(filter(ploidyTable,type=="loss")$sum)
		totalLossSeg = sum(filter(ploidyTable,type=="loss")$`n()`)
		
		
		#SUMMARY STATS
		row <- data.frame(ID=file,ploidy=ploidy,gainSeg=totalGainSeg,totalGain=totalGain,lossSeg=totalLossSeg,totalLoss=totalLoss,file=file)
		cnvPerCellSummary <- rbind(row,cnvPerCellSummary)
		
		#SUMMARIZE INDIVIDUAL EVENTS
		df = tmp %>% select(c(seqnames,start,end,width,name,type))
		if (nrow(df)>0){
			df$file=file
			df$ploidy=ploidy
			df$gene=file
			cnvPerCell<- rbind(df,cnvPerCell)
		}
	}
	
	write.table(cnvPerCell,file = paste0("/Users/zeidh/Desktop/THESIS/CNAs/",cutoff/10000,"cnvPerCell.txt"),quote = F,row.names = F,col.names = T,sep = "\t")
	write.table(cnvPerCellSummary,paste0("/Users/zeidh/Desktop/THESIS/CNAs/",cutoff/10000,"cnvPerCellSummary.txt"),quote = F,row.names = F,col.names = T,sep = "\t")
	
}

copyNumberFinder(10000000,file)



######################################################################
####   Duplications   ####
######################################################################

cnvPerCell=read.table("/Users/zeidh/Desktop/THESIS/CNAs/1000cnvPerCell.txt",header=T)
gains <- dplyr::filter(cnvPerCell,cnvPerCell$type=="gain")
gains$file<-as.factor(gains$file)


cnvpath=paste0("/Users/zeidh/Desktop/THESIS/CNAs/Duplications/")
if (!file.exists(cnvpath) ) { dir.create(cnvpath)}

chr=levels(droplevels(gains$seqnames))[1]
for (chr in levels(droplevels(gains$seqnames)) ){
	filepath=paste0(cnvpath,"/",chr,"/")
	if (!file.exists(filepath) ) { dir.create(filepath)}
	
	breaks = filter(gains,seqnames==chr)
	if (nrow(breaks)>50){
		#breakpointR::plotBreakpoints(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data_pe//",levels(droplevels(breaks$file)),".RData"))
		plotBreakpointsPerChr2.0(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data_pe//",levels(droplevels(breaks$file)),".RData")[1:50], plotspath = filepath,chromosomes = chr)
		}
		
	else 
		{plotBreakpointsPerChr2.0(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data_pe//",levels(droplevels(breaks$file)),".RData"), plotspath = filepath,chromosomes = chr)}
}

for (chr in levels(droplevels(gains$seqnames)) ){
	filepath=paste0(cnvpath,"/",chr,"/")
	if (!file.exists(filepath) ) { dir.create(filepath)}
	
	breaks = dplyr::filter(gains,seqnames==chr)
	gene = read.table("/Users/zeidh/Desktop/THESIS/05-pe.se.metrics.good.quality.ploidy.txt",header=T)%>% select(c(file,gene))
	breaks=merge(breaks,gene,by="file")
	
	write.table(breaks,paste0(filepath,"CNAs.txt"),quote = F,col.names = T,row.names = F,sep = "\t")
}

######################################################################
####   Deletions  ####
######################################################################


deletions <- dplyr::filter(cnvPerCell,cnvPerCell$type=="loss")
deletions$file<-as.factor(deletions$file)

cnvpath=paste0("/Users/zeidh/Desktop/THESIS/CNAs/Deletions/")
if (!file.exists(cnvpath) ) { dir.create(cnvpath)}


for (chr in levels(droplevels(deletions$seqnames)) ){
	filepath=paste0(cnvpath,"/",chr,"/")
	if (!file.exists(filepath) ) { dir.create(filepath)}
	
	breaks = dplyr::filter(deletions,seqnames==chr)
	
	if (nrow(breaks)>100){
		#breakpointR::plotBreakpoints(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data_pe//",levels(droplevels(breaks$file)),".RData"))
		plotBreakpointsPerChr2.0(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data_pe//",levels(droplevels(breaks$file)),".RData")[1:100], plotspath = filepath,chromosomes = chr)
	}
	
	else 
	{plotBreakpointsPerChr2.0(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data_pe//",levels(droplevels(breaks$file)),".RData"), plotspath = filepath,chromosomes = chr)}
}


for (chr in levels(droplevels(deletions$seqnames)) ){
	filepath=paste0(cnvpath,"/",chr,"/")
	if (!file.exists(filepath) ) { dir.create(filepath)}
	
	breaks = dplyr::filter(deletions,seqnames==chr)
	gene = read.table("/Users/zeidh/Desktop/THESIS/05-pe.se.metrics.good.quality.ploidy.txt",header=T)%>% select(c(file,gene))
	breaks=merge(breaks,gene,by="file")
	
	write.table(breaks,paste0(filepath,"CNAs.txt"),quote = F,col.names = T,row.names = F,sep = "\t")
}


