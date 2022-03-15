#################################################################################
#################################################################################
library(GenomicRanges)
library(tidyverse)

#################################################################################
#################################################################################

essentialGenesPlus=read.table("Enrichment/Intersecting_BED/plusEssentialGenes_KBM.bed",header=F)
colnames(essentialGenesPlus)=c("seqnames", "start" ,"end", "V4", "V5")
essentialGenesPlus=essentialGenesPlus[essentialGenesPlus$seqnames %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),]
essentialGenesPlus$seqnames=paste0("chr",essentialGenesPlus$seqnames)
gr_ess=GRanges(essentialGenesPlus)

nonEssentialGenesPlus=read.table("Enrichment/Intersecting_BED/plusNonEssentialGenes_KBM.bed",header=F)
colnames(nonEssentialGenesPlus)=c("seqnames", "start" ,"end", "V4", "V5")
nonEssentialGenesPlus=nonEssentialGenesPlus[nonEssentialGenesPlus$seqnames %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),]
nonEssentialGenesPlus$seqnames=paste0("chr",nonEssentialGenesPlus$seqnames)
gr_non_ss=GRanges(nonEssentialGenesPlus)

G4s=read.table("Enrichment/Intersecting_BED/G4_K.bed")
colnames(G4s)=c("seqnames", "start" ,"end", "V4")
gr_g4=GRanges(G4s)

m=read.table("Enrichment/TRC/CD_HO.bed",header=F)
colnames(m)=c( "seqnames" ,"start" ,"end" ,"V4")
m=m[m$seqnames %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),]
m$seqnames=paste0("chr",m$seqnames)
CD_HO=GRanges(m)

m=read.table("Enrichment/TRC/CD_CD.bed",header=F)
colnames(m)=c( "seqnames" ,"start" ,"end" ,"V4")
m=m[m$seqnames %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),]
m$seqnames=paste0("chr",m$seqnames)
CD_CD=GRanges(m)

m=read.table("Enrichment/TRC/HO_HO.bed",header=F)
colnames(m)=c( "seqnames" ,"start" ,"end" ,"V4")
m=m[m$seqnames %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),]
m$seqnames=paste0("chr",m$seqnames)
HO_HO=GRanges(m)

m=read.table("Enrichment/TRC/HO_CD.bed",header=F)
colnames(m)=c( "seqnames" ,"start" ,"end" ,"V4")
m=m[m$seqnames %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),]
m$seqnames=paste0("chr",m$seqnames)
HO_CD=GRanges(m)


#################################################################################
#################################################################################

calcPercGenesOverlap = function(maxgap,gr1,gr2,ori){
	ess=suppressWarnings(  (length(gr1[queryHits(findOverlaps(gr1,ori, type="any",maxgap = maxgap)),]) )/length(gr_ess) )* 100
	nes=suppressWarnings(  (length(gr2[queryHits(findOverlaps(gr2,ori, type="any",maxgap = maxgap)),]) )/length(gr_non_ss) )* 100
	message("There are ",round(ess,digits = 2),"% of essential genes overlapping with ",deparse(substitute(ori)))
	message("There are ",round(nes,digits = 2),"% of non-essential genes overlapping with ",deparse(substitute(ori)))
}
calcPercGenesOverlap(100000,gr_ess,gr_non_ss,CD_HO)
calcPercGenesOverlap(100000,gr_ess,gr_non_ss,CD_CD)
calcPercGenesOverlap(100000,gr_ess,gr_non_ss,HO_HO)
calcPercGenesOverlap(100000,gr_ess,gr_non_ss,HO_CD)

calcPercGenesOverlap(10000,gr_ess,gr_non_ss,CD_HO)
calcPercGenesOverlap(10000,gr_ess,gr_non_ss,CD_CD)
calcPercGenesOverlap(10000,gr_ess,gr_non_ss,HO_HO)
calcPercGenesOverlap(10000,gr_ess,gr_non_ss,HO_CD)

calcPercGenesOverlap(0,gr_ess,gr_non_ss,CD_HO)
calcPercGenesOverlap(0,gr_ess,gr_non_ss,CD_CD)
calcPercGenesOverlap(0,gr_ess,gr_non_ss,HO_HO)
calcPercGenesOverlap(0,gr_ess,gr_non_ss,HO_CD)

calcPercGenesLeftRight = function(gr,ori,gap){
	left=0
	right=0
	not_intersecting=0
	for (el in 1:length(gr)){
		tmp=gr[el]
		if (length( suppressWarnings( tmp[queryHits(findOverlaps(tmp,ori, type="any",maxgap = gap)),])) > 0) {
			origin=suppressWarnings( ori[queryHits(findOverlaps(ori,tmp, type="any",maxgap = gap)),])
			if (length(origin)<2){
				if (origin@ranges@start > tmp@ranges@start){
					left=left+1} 
				else{right =right+1}
			}else{
				for (i in 1:length(origin)){
					tmp2=origin[i]
					if (tmp2@ranges@start > tmp@ranges@start){left=left+1} 
					else{right =right+1}
			}}
		} else{ not_intersecting=not_intersecting+1}
	}
	message("Check #1: ",left+right+not_intersecting," = ",length(gr))
	message("Check #2: ",left+right," = ",suppressWarnings( sum(countOverlaps(gr,ori, type="any",maxgap = gap))))
	message("Which is ",round(((left+right)/(left+right+not_intersecting))*100,digits = 2),"% ")
	message("That is ",left, " on the left and ", right," on the right for ",deparse(substitute(ori)))
}

#essential genes
calcPercGenesLeftRight(gr=gr_ess,ori=CD_CD,gap=10000)
calcPercGenesLeftRight(gr=gr_ess,ori=CD_HO,gap=10000)
calcPercGenesLeftRight(gr=gr_ess,ori=HO_HO,gap=10000)
calcPercGenesLeftRight(gr=gr_ess,ori=HO_CD,gap=10000)

# non essential genes
calcPercGenesLeftRight(gr=gr_non_ss,ori=CD_CD,gap=10000)
calcPercGenesLeftRight(gr=gr_non_ss,ori=CD_HO,gap=10000)
calcPercGenesLeftRight(gr=gr_non_ss,ori=HO_HO,gap=10000)
calcPercGenesLeftRight(gr=gr_non_ss,ori=HO_CD,gap=10000)


