

BiocManager::install("regioneR")
library(regioneR)
library(tidyverse)
library(GenomicRanges)
library(dplyr)


########################################################################################
################################ Read and save SCE bed files  ########################
########################################################################################

levels=c()
for (file in list.files("/Users/zeidh/Desktop/THESIS/Enrichment/SCEs/ByGene/",full.names = T)){
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

centromeres=read.table("/Users/zeidh/Desktop/THESIS/Anxilliary/centromeres2.txt",header=F)
colnames(centromeres)=c("chr","start","end")
mask=GRanges(centromeres)


########################################################################################
################################ Permutation enrichment analysis #######################
########################################################################################
essentialGenes=read.table("Enrichment/Genes/Initial/essentialKBM7genes.txt",header = T)
essential=merge(genes,essentialGenes,by.y="Gene",by.x="geneName") %>% dplyr::select(chrom,chromStart,chromEnd,strand,transcriptType,geneName)
plusGenes = filter(essential,strand=="+")
minusGenes = filter(essential,strand=="-")
plusProteinGenes = filter(plusGenes,transcriptType=="protein_coding")
gr_plusProteinGenes=GRanges(plusProteinGenes)
minusProteinGenes = filter(minusGenes,transcriptType=="protein_coding")
plus_genes=read.table("Enrichment/Genes/plusEssentialGenes_KBM.bed",header=F) 
colnames(plus_genes)=c("seqnames","start","end","v4","v5")

plus_genes_with_g4=read.table("Enrichment/Intersecting_BED/essentialGenes_KBM_Plus_G4_K.bed",header=F)
colnames(plus_genes_with_g4)=c("seqnames","start","end","v4","v5")
plus_genes_without_g4=read.table("Enrichment/Intersecting_BED/essentialGenes_KBM_without_Plus_G4_K.bed",header=F)

plus_genes=GRanges(plus_genes)

G4s=read.table("Enrichment/Intersecting_BED/G4_K.bed")
colnames(G4s)=c("seqnames", "start" ,"end", "V4")
gr_g4=GRanges(G4s)

plus_genes_with_g4=suppressWarnings(subsetByOverlaps(plus_genes,gr_g4))
plus_genes_without_g4=suppressWarnings(subsetByOverlaps(plus_genes,gr_g4,invert = T))
plusProteinGenes_g4=suppressWarnings(subsetByOverlaps(gr_plusProteinGenes,gr_g4))


nrow(plusGenes)
nrow(plusProteinGenes)
nrow(plus_genes)
length(plus_genes)
nrow(plus_genes_with_g4)
length(plus_genes_with_g4)
nrow(plus_genes_without_g4)
length(plus_genes_without_g4)

B=as.data.frame(plusProteinGenes_g4)
colnames(B)=c("seqnames","start","end","v3","v4","v5","v6")
B=B[B$seqnames %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),]
B$seqnames=paste0("chr",B$seqnames)
B=GRanges(B)

null=data.frame()
true=data.frame()
pval=data.frame()

for (level in levels){
	print(level)
	if (level == "BLM"){
		tmp = sample_n(get(level),size = 1500)
	} else {tmp=get(level)}
	A=GRanges(tmp)
	set.seed(123)
	pt <- overlapPermTest(A=A, B=B, ntimes=100,mc.set.seed=FALSE,force.parallel=TRUE,mask=mask,count.once=TRUE)
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
null=null %>% mutate(color = ifelse(gene == "WT" |gene == "NA878", yes = "#89E088",no= "#88BCE0"))


null$gene=factor(null$gene,levels = c("WT","WRN_RECQL5",  "WRN"  , "RECQL5"  ,"BLM"    ,   "BLM_RECQL5"   ,"NA878"     ) )

ggplot(null)+geom_violin(aes(gene,enrich,fill=color),lwd=0.9,trim=F) +
	geom_boxplot(width=0.05,aes(gene,enrich),lwd=0.9) +
	theme_classic(base_size = 15) +
	geom_point(data=true,aes(gene,enrich),fill="red",colour="black",pch=21, size=5)+
	theme(legend.position = "none",aspect.ratio = 0.41)  +
	labs(x="Cells",y="Enrichment",title = "(+)Essential Genes with G4")+
	geom_text(data=pval,aes(x=gene,y=(max(true$enrich,null$enrich)+0.2),label=sig))

