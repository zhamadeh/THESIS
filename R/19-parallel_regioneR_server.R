########################################################################################
########################################### Packages #################################
########################################################################################

args=commandArgs(trailingOnly = T)
#args[1]=list.files("ALL/",full.names=T)[1]
BiocManager::install("regioneR")
library(regioneR)
library(tidyverse)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
options(bitmapType="cairo")


########################################################################################
################################ Read and save SCE bed files  ########################
########################################################################################

levels=c()
for (file in list.files("SCEs/",full.names = T)){
	print(file)
	tmp=read.table(file,header=F)
	tmp$V4=tmp$V3-tmp$V2
	tmp= filter(tmp,V4<25000)
	colnames(tmp)=c("chr","start","end","width")
	tmp$chr=paste0("chr",tmp$chr)
	assign(strsplit(x=basename(file),split = "[.]")[[1]][1],tmp)
	levels=append(x = levels,values = strsplit(x=basename(file),split = "[.]")[[1]][1])
}

########################################################################################
################################ Centromeres as masking region ########################
########################################################################################

centromeres=read.table("Anxilliary/centromeres2.txt",header=F)
colnames(centromeres)=c("chr","start","end")
mask=GRanges(centromeres)


########################################################################################
################################ Permutation enrichment analysis #######################
########################################################################################

message("Starting on ",basename(args[1])," ... " )

B=read.table(args[1])
colnames(B)[1:3]=c("seqnames","start","end")
B=B[B$seqnames %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),]
B$seqnames=paste0("chr",B$seqnames)
B=GRanges(B)

null=data.frame()
true=data.frame()
pval=data.frame()

for (level in levels){
	print(level)
	#if (level == "BLM"){
	#	tmp = sample_n(get(level),size = 1500)
	#} else {}
	tmp=get(level)
	A=GRanges(tmp)
	set.seed(123)
	pt <- overlapPermTest(A=A, B=B, ntimes=1000,mc.set.seed=FALSE,force.parallel=TRUE,mask=mask)
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
null=null %>% mutate(color = ifelse(gene == "WT" | gene == "NA878", yes = "#89E088",no= "#88BCE0"))

null$gene=factor(null$gene,levels = c("WT","WRN_RECQL5", "RECQL1", "WRN"  , "RECQL5"  ,"BLM"    ,   "BLM_RECQL5"   ,"NA878"     ) )


write.table(null,paste0("Plots/Data/",strsplit(x=basename(args[1]),"[.]")[[1]][1],"_null.txt"),quote = F,row.names = F,col.names = T,sep="\t")
write.table(true,paste0("Plots/Data/",strsplit(x=basename(args[1]),"[.]")[[1]][1],"_true.txt"),quote = F,row.names = F,col.names = T,sep="\t")
write.table(pval,paste0("Plots/Data/",strsplit(x=basename(args[1]),"[.]")[[1]][1],"_pval.txt"),quote = F,row.names = F,col.names = T,sep="\t")



ggplot(null)+geom_violin(aes(gene,enrich,fill=color),lwd=0.9,trim=F) +
	geom_boxplot(width=0.05,aes(gene,enrich),lwd=0.9) +
	theme_classic(base_size = 15) +
	geom_point(data=true,aes(gene,enrich),fill="red",colour="black",pch=21, size=5)+
	theme(legend.position = "none")+#,aspect.ratio = 0.41)+
	labs(x="Cells",y="Enrichment",title = basename)+
	scale_x_discrete(breaks = c("WT","WRN_RECQL5", "RECQL1", "WRN"  , "RECQL5"  ,"BLM"    ,   "BLM_RECQL5"   ,"NA878"   ),labels = c("WT","WRN\nRECQL5", "RECQL1", "WRN"  , "RECQL5"  ,"BLM"    ,   "BLM\nRECQL5"   ,"NA878"   ))+
	geom_text(data=pval,aes(x=gene,y=(max(true$enrich,null$enrich)+0.2),label=sig))  +
	ggsave(paste0("Plots/",basename,".png"))


