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
		
		
		## BUILD PLOIDY TABLE SUMMARY TO CLASSIFY NATIVE STATE AND deletions & LOSSES
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
		
		
		## TURN PLOIDY INTO NUMERIC AND CLASSIFY  deletions/LOSSES
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
	
	write.table(cnvPerCell,file = "/Users/zeidh/Desktop/THESIS/CNAs/20Mb_cnvPerCell.txt",quote = F,row.names = F,col.names = T,sep = "\t")
	write.table(cnvPerCellSummary,file="/Users/zeidh/Desktop/THESIS/CNAs/20Mb_cnvPerCellSummary.txt",quote = F,row.names = F,col.names = T,sep = "\t")
	
}

file = "/Users/zeidh/Desktop/THESIS/CNAs/Aneufinder/binsize_1e+06_PE_CNV.bed.gz"
copyNumberFinder(20000000,file)

file = "/Users/zeidh/Desktop/THESIS/CNAs/Aneufinder/binsize_1e+06_SE_CNV.bed.gz"
copyNumberFinder(20000000,file)



######################################################################
####   Duplications   ####
######################################################################

cnvPerCell1=read.table("/Users/zeidh/Desktop/THESIS/CNAs/Data/01-20Mb_pe_cnvPerCell.txt",header=T)
cnvPerCell2=read.table("/Users/zeidh/Desktop/THESIS/CNAs/Data/01-20Mb_se_cnvPerCell.txt",header=T)

cnvPerCell=rbind(cnvPerCell1,cnvPerCell2)

gains <- dplyr::filter(cnvPerCell,cnvPerCell$type=="gain")
gains$file<-as.factor(gains$file)

gains$library <- gains$file
suppressWarnings(gains <- gains %>% separate(library, c("a","b","c","d","e","f"), "[_-]+"))
gains$gene <- "gene"
for (row in 1:nrow(gains)){
	
	for (letter in c("a","b","c","d","e","f")){
		if (is.na(gains[row,letter])!=T){
			if (gains[row,letter]=="WT" | gains[row,letter]=="wt"){
				gains[row,"gene"]="WT"
			}
			else if (gains[row,letter]=="blm" | gains[row,letter]=="BLM"  | gains[row,letter]=="blm1" ) {
				gains[row,"gene"]="BLM"
			}
			else if (gains[row,letter]=="wrn"  ) {
				if (gains[row,"a"]=="recql5"){
					gains[row,"gene"]="WRN/RECQL5"
				}
				else{
					gains[row,"gene"]="WRN"
				}
			}
			
			else if (gains[row,letter]=="RECQL5" | gains[row,letter]=="recql5" | gains[row,letter]=="RECQ5" | gains[row,letter]=="recq5" ) {
				if (gains[row,"gene"]=="BLM"){
					gains[row,"gene"]="BLM/RECQL5"
				}
				else{
					gains[row,"gene"]="RECQL5"
				}
			}
			else if (gains[row,letter]=="blmrecq5"){
				gains[row,"gene"]="BLM/RECQL5"
			}				
			else if (gains[row,letter]=="RECQL1"| gains[row,letter]=="recql1"){
				gains[row,"gene"]="RECQL1"
			}
			else if ( gains[row,letter]=="fucciwt"| gains[row,letter]=="kbm7" | gains[row,letter]=="wtfucci" ){
				gains[row,"gene"]="WT"
			}
			else if (gains[row,letter]=="RTEL"| gains[row,letter]=="rtel"){
				gains[row,"gene"]="RTEL1"
			}
		}
	}
	
}
gains <- select(gains,-c(a,b,c,d,e,f))
gains$gene=as.factor(gains$gene)




wt = filter(gains,gene=="WT")
WT=GRanges(wt)

wtpath=paste0("/Users/zeidh/Desktop/THESIS/CNAs/WT/")
if (!file.exists(wtpath) ) { dir.create(wtpath)}
for (chr in levels(droplevels(wt$seqnames)) ){
	wt_chr=filter(wt,seqnames==chr)
	filepath2=paste0(wtpath,"/",chr,"/")
	if (!file.exists(filepath2) ) { dir.create(filepath2)}
	plotBreakpointsPerChr2.0(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data/",levels(droplevels(wt_chr$file)),".RData"), plotspath = filepath2,chromosomes = chr,tmp=wt_chr)
}




gains=gains[gains$gene!="WT",]
cnvpath=paste0("/Users/zeidh/Desktop/THESIS/CNAs/Somatic/Duplications/")
if (!file.exists(cnvpath) ) { dir.create(cnvpath)}
fucci_breaks_dup = data.frame()

for (chr in levels(droplevels(gains$seqnames)) ){
	filepath=paste0(cnvpath,"/",chr,"/")
	filepath2=paste0(filepath,"/control_ko/")
	filepath3=paste0(filepath,"/control_wt/")
	
	breaks = filter(gains,seqnames==chr)
	wt_chr=filter(wt,seqnames==chr)
	if (nrow(wt_chr)>0){
		WT_chr=GRanges(wt_chr)
	} else{
		next
	}

	overlaps = countOverlaps(GRanges(breaks),WT_chr, type="any",minoverlap = 0.5)
	overlap_rows = which(overlaps > 2)
	overlaps2 = countOverlaps(WT_chr,GRanges(breaks), type="any",minoverlap = 0.5)
	overlap_rows2 = which(overlaps2 > 2)
	tmp = breaks[-overlap_rows,]
	tmp2 = breaks[overlap_rows,]
	tmp3 = wt_chr[overlap_rows2,]
	#tmp4=rbind(tmp2,tmp3)
	
	if (nrow(tmp)>0){
		print(chr)
		fucci_breaks_dup<- rbind(fucci_breaks_dup,tmp)
		if (!file.exists(filepath) ) { dir.create(filepath)}
		plotBreakpointsPerChr2.0(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data/",levels(droplevels(tmp$file)),".RData"), plotspath = filepath,chromosomes = chr,tmp=tmp)
		if (nrow(tmp2)>0){
			if (!file.exists(filepath2) ) { dir.create(filepath2)}
			if (nrow(tmp2)<10){
				plotBreakpointsPerChr2.0(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data/",levels(droplevels(tmp2$file)),".RData"), plotspath = filepath2,chromosomes = chr,tmp=breaks)
			}  else{ plotBreakpointsPerChr2.0(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data/",sample(levels(droplevels(tmp2$file)),size=10),".RData"), plotspath = filepath2,chromosomes = chr,tmp=breaks)
			}
		}
		if (nrow(tmp3)>0){
			if (!file.exists(filepath3) ) { dir.create(filepath3)}
			if (nrow(tmp3)<10){
				plotBreakpointsPerChr2.0(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data/",levels(droplevels(tmp3$file)),".RData"), plotspath = filepath3,chromosomes = chr,tmp=wt_chr)
			}  else{ plotBreakpointsPerChr2.0(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data/",sample(levels(droplevels(tmp3$file)),size=10),".RData"), plotspath = filepath3,chromosomes = chr,tmp=wt_chr)
			}
		}
	} 
}

write.table(fucci_breaks_dup,"CNAs/Data/04-somatic_duplications.txt",quote = F,row.names = F,col.names = T,sep = "\t")

######################################################################
####   Deletions  ####
######################################################################


deletions <- dplyr::filter(cnvPerCell,cnvPerCell$type=="loss")
deletions$file<-as.factor(deletions$file)

deletions$library <- deletions$file
suppressWarnings(deletions <- deletions %>% separate(library, c("a","b","c","d","e","f"), "[_-]+"))
deletions$gene <- "gene"
for (row in 1:nrow(deletions)){
	
	for (letter in c("a","b","c","d","e","f")){
		if (is.na(deletions[row,letter])!=T){
			if (deletions[row,letter]=="WT" | deletions[row,letter]=="wt"){
				deletions[row,"gene"]="WT"
			}
			else if (deletions[row,letter]=="blm" | deletions[row,letter]=="BLM"  | deletions[row,letter]=="blm1" ) {
				deletions[row,"gene"]="BLM"
			}
			else if (deletions[row,letter]=="wrn"  ) {
				if (deletions[row,"a"]=="recql5"){
					deletions[row,"gene"]="WRN/RECQL5"
				}
				else{
					deletions[row,"gene"]="WRN"
				}
			}
			
			else if (deletions[row,letter]=="RECQL5" | deletions[row,letter]=="recql5" | deletions[row,letter]=="RECQ5" | deletions[row,letter]=="recq5" ) {
				if (deletions[row,"gene"]=="BLM"){
					deletions[row,"gene"]="BLM/RECQL5"
				}
				else{
					deletions[row,"gene"]="RECQL5"
				}
			}
			else if (deletions[row,letter]=="blmrecq5"){
				deletions[row,"gene"]="BLM/RECQL5"
			}				
			else if (deletions[row,letter]=="RECQL1"| deletions[row,letter]=="recql1"){
				deletions[row,"gene"]="RECQL1"
			}
			else if ( deletions[row,letter]=="fucciwt"| deletions[row,letter]=="kbm7" | deletions[row,letter]=="wtfucci" ){
				deletions[row,"gene"]="WT"
			}
			else if (deletions[row,letter]=="RTEL"| deletions[row,letter]=="rtel"){
				deletions[row,"gene"]="RTEL1"
			}
		}
	}
	
}
deletions <- select(deletions,-c(a,b,c,d,e,f))
deletions$gene=as.factor(deletions$gene)




wt = filter(deletions,gene=="WT")
WT=GRanges(wt)

wtpath=paste0("/Users/zeidh/Desktop/THESIS/CNAs/WT/")
if (!file.exists(wtpath) ) { dir.create(wtpath)}
for (chr in levels(droplevels(wt$seqnames)) ){
	wt_chr=filter(wt,seqnames==chr)
	filepath2=paste0(wtpath,"/",chr,"/")
	if (!file.exists(filepath2) ) { dir.create(filepath2)}
	plotBreakpointsPerChr2.0(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data/",levels(droplevels(wt_chr$file)),".RData"), plotspath = filepath2,chromosomes = chr,tmp=wt_chr)
}




deletions=deletions[deletions$gene!="WT",]
cnvpath=paste0("/Users/zeidh/Desktop/THESIS/CNAs/Somatic/Deletions/")
if (!file.exists(cnvpath) ) { dir.create(cnvpath)}
fucci_breaks_deletions = data.frame()

for (chr in levels(droplevels(deletions$seqnames)) ){
	filepath=paste0(cnvpath,"/",chr,"/")
	filepath2=paste0(filepath,"/control_ko/")
	filepath3=paste0(filepath,"/control_wt/")
	
	breaks = filter(deletions,seqnames==chr)
	wt_chr=filter(wt,seqnames==chr)
	if (nrow(wt_chr)>0){
		WT_chr=GRanges(wt_chr)
	} else{
		next
	}

	overlaps = countOverlaps(GRanges(breaks),WT_chr, type="any",minoverlap = 0.5)
	overlap_rows = which(overlaps > 2)
	overlaps2 = countOverlaps(WT_chr,GRanges(breaks), type="any",minoverlap = 0.5)
	overlap_rows2 = which(overlaps2 > 2)
	tmp = breaks[-overlap_rows,]
	tmp2 = breaks[overlap_rows,]
	tmp3 = wt_chr[overlap_rows2,]
	#tmp4=rbind(tmp2,tmp3)
	
	if (nrow(tmp)>0){
		print(chr)
		fucci_breaks_deletions<- rbind(fucci_breaks_deletions,tmp)
		if (!file.exists(filepath) ) { dir.create(filepath)}
		#plotBreakpointsPerChr2.0(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data/",levels(droplevels(tmp$file)),".RData"), plotspath = filepath,chromosomes = chr,tmp=tmp)
		if (nrow(tmp2)>0){
			if (!file.exists(filepath2) ) { dir.create(filepath2)}
			if (nrow(tmp2)<10){
				#plotBreakpointsPerChr2.0(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data/",levels(droplevels(tmp2$file)),".RData"), plotspath = filepath2,chromosomes = chr,tmp=breaks)
			}  else{ #plotBreakpointsPerChr2.0(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data/",sample(levels(droplevels(tmp2$file)),size=10),".RData"), plotspath = filepath2,chromosomes = chr,tmp=breaks)
			}
		}
		if (nrow(tmp3)>0){
			if (!file.exists(filepath3) ) { dir.create(filepath3)}
			if (nrow(tmp3)<10){
				#plotBreakpointsPerChr2.0(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data/",levels(droplevels(tmp3$file)),".RData"), plotspath = filepath3,chromosomes = chr,tmp=tmp3)
			}  else{ #plotBreakpointsPerChr2.0(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data/",sample(levels(droplevels(tmp3$file)),size=10),".RData"), plotspath = filepath3,chromosomes = chr,tmp=wt_chr)
			}
		}
	} 
}

write.table(fucci_breaks_deletions,"CNAs/Data/03-somatic_deletions.txt",quote = F,row.names = F,col.names = T,sep = "\t")

######################################################################
####   Somatic deletions and duplications ####
######################################################################

write.table(rbind(fucci_breaks_dup,fucci_breaks_deletions),"CNAs/Data/02-somatic_cnas.txt",quote = F,row.names = F,col.names = T,sep = "\t")


######################################################################
####   WT Deletions ####
######################################################################

dels<- dplyr::filter(cnvPerCell,cnvPerCell$type=="loss")
wt = filter(dels,gene=="WT")
WT=GRanges(wt)

somatic_rows=c()
for (i in 1:length(WT)){
	gr1= WT[i]
	gr2=  WT[-i]
	
	overlaps = countOverlaps(gr1,gr2, type="any",minoverlap = 0.5)
	if (overlaps < 3){
		somatic_rows=append(i,somatic_rows)
	}
}
wt_somatic= wt[somatic_rows,]
wtpath=paste0("/Users/zeidh/Desktop/THESIS/CNAs/WT/Somatic/")
if (!file.exists(wtpath) ) { dir.create(wtpath)}
for (chr in levels(droplevels(wt_somatic$seqnames))){
	wt_chr=filter(wt_somatic,seqnames==chr)
	filepath2=paste0(wtpath,"/",chr,"/")
	if (nrow(wt_chr)>0){
		if (!file.exists(filepath2) ) { dir.create(filepath2)}
		plotBreakpointsPerChr2.0(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data/",levels(droplevels(wt_chr$file)),".RData"), plotspath = filepath2,chromosomes = chr,tmp=wt_chr)
		
	}
}
a=wt_somatic[wt_somatic$seqnames=="chr1" & wt_somatic$file=="wtfucci-hyd-a1h5-f10-r1-c5_S178_.trimmed.mdup.bam",]
b=wt_somatic[wt_somatic$seqnames=="chr2" & wt_somatic$file!="kbm7-a1h5-hkt-r40-c22_S130_.trimmed.mdup.bam",]
c=wt_somatic[wt_somatic$seqnames=="chr3" & wt_somatic$file%in% c("kbm7-a1h5-hkt-r40-c34_S142_.trimmed.mdup.bam",
																 "wtfucci-hyd-a1h5-f10-r1-c5_S178_.trimmed.mdup.bam",
																 "wtfucci-hyd-a1h5-f10-r5-c5_S192_.trimmed.mdup.bam"),]
d=wt_somatic[wt_somatic$seqnames=="chr4" & wt_somatic$file!="wt-fucci-fdu-single-f9-r2-c2_S201_.trimmed.mdup.bam",]
e=wt_somatic[wt_somatic$seqnames=="chr9" ,]
f=wt_somatic[wt_somatic$seqname=="chr11" & wt_somatic$file=="wtfucci-hyd-mn-f03-r5-c3_S66_.trimmed.mdup.bam",]
g=wt_somatic[wt_somatic$seqnames=="chr13" & wt_somatic$file=="kbm7-yeast-a1h5-hkt-r16-c67_S445_.trimmed.mdup.bam",]
h=wt_somatic[wt_somatic$seqnames=="chr16" & wt_somatic$file=="fucci-wt-single-f8-r7-c6_S97_.trimmed.mdup.bam",]
i=wt_somatic[wt_somatic$seqnames=="chr20",]

deletions=rbind(a,b,c,d,e,f,g,h,i)
write.table(deletions,"CNAs/Data/03-somatic_deletions_wt.txt",quote = F,row.names = F,col.names = T,sep = "\t")


######################################################################
####   WT ####
######################################################################
gains <- dplyr::filter(cnvPerCell,cnvPerCell$type=="gain")
wt = filter(gains,gene=="WT")
WT=GRanges(wt)

somatic_rows=c()
for (i in 1:length(WT)){
	gr1= WT[i]
	gr2=  WT[-i]
	
	overlaps = countOverlaps(gr1,gr2, type="any",minoverlap = 0.5)
	if (overlaps < 3){
		somatic_rows=append(i,somatic_rows)
	}
}
wt_somatic= wt[somatic_rows,]
wtpath=paste0("/Users/zeidh/Desktop/THESIS/CNAs/WT/Somatic/Duplications/")
if (!file.exists(wtpath) ) { dir.create(wtpath)}
for (chr in levels(droplevels(wt_somatic$seqnames))){
	wt_chr=filter(wt_somatic,seqnames==chr)
	filepath2=paste0(wtpath,"/",chr,"/")
	if (nrow(wt_chr)>0){
		if (!file.exists(filepath2) ) { dir.create(filepath2)}
		plotBreakpointsPerChr2.0(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data/",levels(droplevels(wt_chr$file)),".RData"), plotspath = filepath2,chromosomes = chr,tmp=wt_chr)
		
	}
}
a=wt_somatic[wt_somatic$seqnames=="chr2" & wt_somatic$file=="kbm7-a1h5-hkt-r69-c30_S1182_.trimmed.mdup.bam",]
b=wt_somatic[wt_somatic$seqnames=="chr3" & wt_somatic$file!="wt-fucci-fdu-single-f10-r4-c3_S237_.trimmed.mdup.bam",]
c=wt_somatic[wt_somatic$seqnames=="chr4" & wt_somatic$file!= "wt-fucci-fdu-single-f10-r4-c3_S237_.trimmed.mdup.bam",]
d=wt_somatic[wt_somatic$seqnames=="chr5" & wt_somatic$file=="fucci-wt-single-f7-r4-c3_S24_.trimmed.mdup.bam",]
e=wt_somatic[wt_somatic$seqnames=="chr9" & wt_somatic$file== "fucci-wt-single-f11-r2-c5_S110_.trimmed.mdup.bam",]
f=wt_somatic[wt_somatic$seqname=="chr13" ,]
g=wt_somatic[wt_somatic$seqnames=="chr20",]
h=wt_somatic[wt_somatic$seqnames=="chrX" & wt_somatic$file!="wt-fucci-fdu-single-f10-r4-c3_S237_.trimmed.mdup.bam",]

duplications=rbind(a,b,c,d,e,f,g,h)
write.table(duplications,"CNAs/Data/04-somatic_duplications_wt.txt",quote = F,row.names = F,col.names = T,sep = "\t")

somatic_WT=rbind(duplications,deletions)
write.table(somatic_WT,"CNAs/Data/02-somatic_cnas_wt.txt",quote = F,row.names = F,col.names = T,sep = "\t")


