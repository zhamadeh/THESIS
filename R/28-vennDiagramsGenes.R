################################################################################################

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111272

BiocManager::install("BioVenn")
require(BioVenn)

#https://www.biovenn.nl/index.php

################################################################################################

hap1_genes=read.table("Enrichment/GSE111272_KBM7_repFPKM_table.txt",header=T)
library(mygene)
gene.list = hap1_genes$gene_id
geneNames=getGenes(gene.list, fields='symbol')
geneName=as.data.frame(geneNames)

merge=merge(hap1_genes,geneName,by.x="gene_id",by.y="query")
merge= dplyr::select(merge,c(symbol, ctrl_1, ctrl_2, ctrl_3))
merge=na.omit(merge)
merge=merge[merge$symbol!="Y_RNA",]


merge=dplyr::select(merge,c(symbol, ctrl_1,ctrl_2,ctrl_3 ))
merge$mean=merge$ctrl_1+merge$ctrl_2+merge$ctrl_3
merge$mean=merge$mean/3

sort=merge[order(-merge$mean),] 
low=sort %>% filter(sort$mean==0)
midHigh=sort %>% filter(sort$mean>0)
midHigh=midHigh[order(-midHigh$mean),]
high=midHigh[1:(nrow(midHigh)/2),]
mid=midHigh[(nrow(midHigh)/2):nrow(midHigh),]
reallyHigh=filter(sort,sort$mean>500)



################################################################################################
merge1=dplyr::select(merge,c(symbol, ctrl_1))
merge2=dplyr::select(merge,c(symbol, ctrl_2 ))
merge3=dplyr::select(merge,c(symbol, ctrl_3 ))

sort1=merge1[order(-merge1$ctrl_1),] 
low1=sort1 %>% filter(sort1$ctrl_1==0)
midHigh=sort1 %>% filter(sort1$ctrl_1>0)
midHigh=midHigh[order(-midHigh$ctrl_1),]
high1=midHigh[1:(nrow(midHigh)/2),]
mid1=midHigh[(nrow(midHigh)/2):nrow(midHigh),]
reallyHigh1=filter(sort1,sort1$ctrl_1>500)

sort2=merge2[order(-merge2$ctrl_2),] 
low2=sort2 %>% filter(sort2$ctrl_2==0)
midHigh=sort2 %>% filter(sort2$ctrl_2>0)
midHigh=midHigh[order(-midHigh$ctrl_2),]
high2=midHigh[1:(nrow(midHigh)/2),]
mid2=midHigh[(nrow(midHigh)/2):nrow(midHigh),]
reallyHigh2=filter(sort2,sort2$ctrl_2>500)

sort3=merge3[order(-merge3$ctrl_3),]
low3=sort3 %>% filter(sort3$ctrl_3==0)
midHigh=sort3 %>% filter(sort3$ctrl_3>0)
midHigh=midHigh[order(-midHigh$ctrl_3),]
high3=midHigh[1:(nrow(midHigh)/2),]
mid3=midHigh[(nrow(midHigh)/2):nrow(midHigh),]
reallyHigh3=filter(sort3,sort3$ctrl_3>500)

length(low1$symbol)
length(low2$symbol)
low=data.frame(symbol=intersect(intersect(low1$symbol,low2$symbol),low3$symbol) )
mid=data.frame(symbol=intersect(intersect(mid1$symbol,mid2$symbol),mid3$symbol) )
high=data.frame(symbol=intersect(intersect(high1$symbol,high2$symbol),high3$symbol) )
reallyHigh=data.frame(symbol=intersect(intersect(reallyHigh1$symbol,reallyHigh2$symbol),reallyHigh3$symbol) )


########################################################################################################################
########################################################################################################################


essentialGenes=read.table("Enrichment/Genes/Initial/essentialKBM7genes.txt",header = T)
essential=merge(genes,essentialGenes,by.y="Gene",by.x="geneName") %>% dplyr::select(chrom,chromStart,chromEnd,strand,transcriptType,geneName)
nonessentialGenes=read.table("Enrichment/Genes/Initial/nonEssentialKBM7genes.txt",header = T)
nonessential=merge(genes,nonessentialGenes,by.y="Gene",by.x="geneName") %>% dplyr::select(chrom,chromStart,chromEnd,strand,transcriptType,geneName)


numHighlyTranscribedGenes=levels(droplevels(as.factor(high$symbol)))
numMidlyTranscribedGenes=levels(droplevels(as.factor(mid$symbol)))
numLowlyTranscribedGenes=levels(droplevels(as.factor(low$symbol)))
numReallyHighlyTranscribedGenes=levels(droplevels(as.factor(reallyHigh$symbol)))
numEssentialGenes=levels(droplevels(as.factor(essential$geneName)))
numNonEssentialGenes=levels(droplevels(as.factor(nonessential$geneName)))

r_high_ess=(length(intersect(numReallyHighlyTranscribedGenes,numEssentialGenes)))
high_ess=(length(intersect(numHighlyTranscribedGenes,numEssentialGenes)))
low_ess=(length(intersect(numLowlyTranscribedGenes,numEssentialGenes)))

r_high_n_ess=(length(intersect(numReallyHighlyTranscribedGenes,numNonEssentialGenes)))
high_n_ess=(length(intersect(numHighlyTranscribedGenes,numNonEssentialGenes)))
low_n_ess=(length(intersect(numLowlyTranscribedGenes,numNonEssentialGenes)))


length(numHighlyTranscribedGenes)
length(numLowlyTranscribedGenes)
length(numEssentialGenes)
length(numNonEssentialGenes)

high_n_ess
low_n_ess
high_ess
low_ess
r_high_n_ess

write.table(numReallyHighlyTranscribedGenes,"RNA/ActivelyTranscribedGenes/r_high.txt",quote = F,row.names = F,col.names = F,sep="\t")
write.table(numHighlyTranscribedGenes,"RNA/ActivelyTranscribedGenes/high.txt",quote = F,row.names = F,col.names = F,sep="\t")
write.table(numMidlyTranscribedGenes,"RNA/ActivelyTranscribedGenes/mid.txt",quote = F,row.names = F,col.names = F,sep="\t")
write.table(numLowlyTranscribedGenes,"RNA/ActivelyTranscribedGenes/low.txt",quote = F,row.names = F,col.names = F,sep="\t")
write.table(numEssentialGenes,"RNA/EssentialGenes//essential.txt",quote = F,row.names = F,col.names = F,sep="\t")
write.table(numNonEssentialGenes,"RNA/EssentialGenes/nonessential.txt",quote = F,row.names = F,col.names = F,sep="\t")


biovenn <- draw.venn( numHighlyTranscribedGenes,numMidlyTranscribedGenes,numLowlyTranscribedGenes,xtitle = "Highly\nExpressed Genes",
					  ytitle = "Midly\nExpressed\nGenes", ztitle = "Lowly\nGenes",
					  xt_s=0.5,yt_s=0.5,zt_s=0.5, subtitle="Example diagram",nr_fb=0,nr_s=0)



## ANother RNA dataset
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71284
#rna = read.table("../GSE71284_control_WT2_WT3_C1_C2_vs_truncation_3kb1_3kb2_100kb1_100kb2_CUFFDIFF_gene_exp.diff",header = T)


########################################################################################################################
########################################################################################################################


genes=read.table("/Users/zeidh/Desktop/THESIS/Enrichment/Genes/Initial/genes.txt",header=T)%>% dplyr::select(c(chrom,chromStart,chromEnd,strand,transcriptType,geneName,geneName2))
genes=dplyr::select(genes,c(chrom,chromStart,chromEnd,strand,transcriptType,geneName,geneName2))
genes$chrom<-gsub(x=genes$chrom,pattern = "chr",replacement = "")
genes=genes[genes$chrom %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),]

reallyHighGenes=merge(reallyHigh,genes,by.x="symbol",by.y="geneName")
highGenes=merge(high,genes,by.x="symbol",by.y="geneName")
midGenes=merge(mid,genes,by.x="symbol",by.y="geneName")
lowGenes=merge(low,genes,by.x="symbol",by.y="geneName")

reallyHighGenes=dplyr::select(reallyHighGenes,c( "chrom" ,"chromStart" , "chromEnd","strand","transcriptType" ,"symbol"))
highGenes=dplyr::select(highGenes,c( "chrom" ,"chromStart" , "chromEnd","strand","transcriptType" ,"symbol"))
midGenes=dplyr::select(midGenes,c( "chrom" ,"chromStart" , "chromEnd","strand","transcriptType" ,"symbol"))
lowGenes=dplyr::select(lowGenes,c( "chrom" ,"chromStart" , "chromEnd","strand","transcriptType" ,"symbol"))

write.table(reallyHighGenes,"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/reallyHighGenes.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(reallyHighGenes,strand=="+"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/plusReallyHighGenes.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(reallyHighGenes,strand=="-"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/minusReallyHighGenes.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(midGenes,"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/midGenes.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(midGenes,strand=="+"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/plusMidGenes.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(midGenes,strand=="-"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/minusMidGenes.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(highGenes,"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/highGenes.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(highGenes,strand=="+"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/plusHighGenes.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(highGenes,strand=="-"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/minusHighGenes.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(lowGenes,"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/lowGenes.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(lowGenes,strand=="+"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/plusLowGenes.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(lowGenes,strand=="-"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/minusLowGenes.bed",quote = F,row.names = F,col.names = F,sep = "\t")


########################################################################################################################
########################################################################################################################

########################################################################################################################
########################################################################################################################



write.table(filter(reallyHighGenes,transcriptType=="CDS"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/reallyHighGenesCDS.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(reallyHighGenes,transcriptType=="exon"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/reallyHighGenesexon.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(reallyHighGenes,transcriptType=="five_prime_utr"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/reallyHighGenesfive_prime_utr.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(reallyHighGenes,transcriptType=="start_codon"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/reallyHighGenesstart_codon.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(reallyHighGenes,transcriptType=="stop_codon"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/reallyHighGenesstop_codon.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(reallyHighGenes,transcriptType=="three_prime_utr"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/reallyHighGenesthree_prime_utr.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(reallyHighGenes,transcriptType=="transcript"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/reallyHighGenestranscript.bed",quote = F,row.names = F,col.names = F,sep = "\t")

write.table(filter(highGenes,transcriptType=="CDS"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/highGenesCDS.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(highGenes,transcriptType=="exon"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/highGenesexon.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(highGenes,transcriptType=="five_prime_utr"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/highGenesfive_prime_utr.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(highGenes,transcriptType=="start_codon"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/highGenesstart_codon.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(highGenes,transcriptType=="stop_codon"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/highGenesstop_codon.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(highGenes,transcriptType=="three_prime_utr"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/highGenesthree_prime_utr.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(highGenes,transcriptType=="transcript"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/highGenestranscript.bed",quote = F,row.names = F,col.names = F,sep = "\t")

write.table(filter(midGenes,transcriptType=="CDS"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/midGenesCDS.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(midGenes,transcriptType=="exon"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/midGenesexon.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(midGenes,transcriptType=="five_prime_utr"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/midGenesfive_prime_utr.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(midGenes,transcriptType=="start_codon"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/midGenesstart_codon.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(midGenes,transcriptType=="stop_codon"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/midGenesstop_codon.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(midGenes,transcriptType=="three_prime_utr"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/midGenesthree_prime_utr.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(midGenes,transcriptType=="transcript"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/midGenestranscript.bed",quote = F,row.names = F,col.names = F,sep = "\t")

write.table(filter(lowGenes,transcriptType=="CDS"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/lowGenesCDS.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(lowGenes,transcriptType=="exon"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/lowGenesexon.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(lowGenes,transcriptType=="five_prime_utr"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/lowGenesfive_prime_utr.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(lowGenes,transcriptType=="start_codon"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/lowGenesstart_codon.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(lowGenes,transcriptType=="stop_codon"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/lowGenesstop_codon.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(lowGenes,transcriptType=="three_prime_utr"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/lowGenesthree_prime_utr.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(filter(lowGenes,transcriptType=="transcript"),"/Users/zeidh/Desktop/THESIS/Enrichment/ActivelyTranscribedGenes/lowGenestranscript.bed",quote = F,row.names = F,col.names = F,sep = "\t")


