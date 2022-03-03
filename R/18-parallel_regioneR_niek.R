########################################################################################
########################################### Packages #################################
########################################################################################

args=commandArgs(trailingOnly = T)
#args=c("Enrichment/Intersecting_BED/GRCh37_G4_1-3.bed")
BiocManager::install("regioneR")
library(regioneR)
library(tidyverse)
library(GenomicRanges)
library(dplyr)


########################################################################################
################################ Read and save SCE bed files  ########################
########################################################################################

levels=c()
for (file in list.files("/Users/zeidh/Desktop/THESIS/SCEs/BED/NIEK/",full.names = T)){
	tmp=read.table(file,header=T)
	assign(strsplit(x=basename(file),split = "[.]")[[1]][1],tmp)
	levels=append(x = levels,values = strsplit(x=basename(file),split = "[.]")[[1]][1])
}

########################################################################################
################################ Centromeres as masking region ########################
########################################################################################

centromeres=read.table("/Users/zeidh/Desktop/THESIS/Anxilliary/centromeres2.txt",header=F)
colnames(centromeres)=c("chr","start","end")
mask=GRanges(centromeres)

########################################################################################
################################ Permutation enrichment analysis #######################
########################################################################################

message("Starting on ",basename(args[1])," ... " )

B=read.table(args[1])
colnames(B)=c("seqnames","start","end")
B=B[B$seqnames %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),]
B$seqnames=paste0("chr",B$seqnames)
B=GRanges(B)

null=data.frame()
true=data.frame()
pval=data.frame()

for (level in levels){
	A=GRanges(filter(niek,Characteristics.cell.line.==level))
	set.seed(123)
	pt <- overlapPermTest(A=A, B=B, ntimes=100,mc.set.seed=FALSE,force.parallel=TRUE,mask=mask)
	tmp1=data.frame(Overlaps=pt$numOverlaps$permuted,gene=level)
	tmp1$enrich = tmp1$Overlaps/mean(tmp1$Overlaps)
	tmp2=data.frame(Overlaps=pt$numOverlaps$observed,gene=level)
	tmp2$enrich=tmp2$Overlaps/mean(tmp1$Overlaps)
	null=rbind(tmp1,null)
	true=rbind(tmp2,true)
	p = data.frame(gene=level,pval=pt$numOverlaps$pval)
	pval=rbind(pval,p)
}

pval=pval %>% mutate(sig = case_when(pval <= 0.001 ~ '***',pval < 0.01 ~ '**',pval < 0.05 ~ '*',pval < 1 ~ 'ns'))

########################################################################################
################################# Plotting #############################################
########################################################################################

null=null %>% mutate(color = ifelse(gene %in% c("WT1","WT2","WT3","WT4") , "#89E088", "#88BCE0"))

null$gene <- factor(null$gene, levels = c( "WT1" ,"WT2" ,"WT3", "WT4","BS1" ,"BS2" ,"BS3", "BS4"))


ggplot(null)+geom_violin(aes(gene,enrich,fill=color),lwd=0.9,trim=F) +
	geom_boxplot(width=0.05,aes(gene,enrich),lwd=0.9) +
	theme_classic(base_size = 15) +
	geom_point(data=true,aes(gene,enrich),fill="red",colour="black",pch=21, size=5)+
	theme(legend.position = "none")  +
	labs(x="Cells",y="Enrichment",title = strsplit(x=basename(args[1]),"[.]")[[1]][1])+
	geom_text(data=pval,aes(x=gene,y=(max(true$enrich,null$enrich)+0.2),label=sig))  +
	ggsave(paste0("Enrichment/EnrichmentPlots/",strsplit(x=basename(args[1]),"[.]")[[1]][1],".png"))


