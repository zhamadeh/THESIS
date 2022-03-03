########################################################################################
########################################### Packages #################################
########################################################################################
# This script estimate the precision recall of BPR in calling SCEs with and without breakpoint processing
# It uses Niek's SCEs from BAIT that he curated as a golden standard benchmark to compare the breakpoints from BPR
# Its a rough estimate since the coordinates from BAIT are likely different than what BPR calculates

library(tidyverse)
library(bedr)
library(GenomicRanges)
library(data.table)

########################################################################################
################################ Niek's bed files and metadata ########################
########################################################################################

metadata <- read.table("Anxilliary/Metadata/E-MTAB-5976.sdrf.txt",sep="\t",header=T)
metadata <-  select(metadata,c(Comment.ENA_RUN., Extract.Name,Factor.Value.protocol., Characteristics.genotype., Characteristics.organism.))

#  Put all bed files in one file
bed <- data.frame()
for (file in list.files("Anxilliary/Benchmark_sce_files/",full.names = T)){
	tmp <- read.table(file,sep="\t")
	if (ncol(tmp)==4){
		tmp$V5=tmp$V3 - tmp$V2
		bed <- rbind(tmp,bed)
	} else {bed <- rbind(tmp,bed)}
}
#merge with metadata
merge <- merge(bed,metadata,by.x="V4",by.y="Extract.Name")

#sort into BLM and WT SCEs
wt <- filter(merge,merge$Characteristics.genotype.=="wild type genotype") %>% select(-c(Characteristics.organism., Characteristics.genotype., Factor.Value.protocol.,)) %>% select(c(V1,V2,V3,V5,V4,Comment.ENA_RUN.))
blm <- filter(merge,merge$Characteristics.genotype.=="homozygous Blm knockout")%>% select(-c(Characteristics.organism., Characteristics.genotype., Factor.Value.protocol.))%>% select(c(V1,V2,V3,V5,V4,Comment.ENA_RUN.))
colnames(blm) = c("chr"   ,  "start"     ,   "end"  , "width"    ,    "Extract.Name","Comment.ENA_RUN." )
blm <- blm %>% select(c(chr,start,end,width,Extract.Name,Comment.ENA_RUN.))
colnames(wt) = c("chr"   ,  "start"     ,   "end"  , "width"    ,    "Extract.Name" ,"Comment.ENA_RUN.")
wt <- wt %>% select(c(chr,start,end,width,Extract.Name,Comment.ENA_RUN.))
bait <- rbind(wt,blm)
bait$window=NA
bait <- rbind(blm,wt)

###############################################################################################################
####################################### unfiltered SCEs from BPR 100reads #####################################
###############################################################################################################

bpr_se_100 = read.table("Precision Recall/Input/breakpoints_100_reads_pe.txt",header = T) # %>% select(-c(CI.start,CI.end))
bpr_pe_100 = read.table("Precision Recall/Input/breakpoints_100_reads_se.txt",header = T) #%>% select(-c(CI.start,CI.end))

bind<- rbind(bpr_pe_100,bpr_se_100)

for (i in 1:5){
	bind$filenames<-tools::file_path_sans_ext(bind$filenames)
}
bind$filenames <- gsub( "_", "", as.character(bind$filenames) )
bind$filenames <- as.factor(bind$filenames)

bpr_merge <- merge(bind,metadata,by.x="filenames",by.y="Comment.ENA_RUN.") %>% select(-c( Factor.Value.protocol., Characteristics.organism.))

bpr_wt <- filter(bpr_merge,bpr_merge$Characteristics.genotype.=="wild type genotype") %>% select(-c(start,end,Characteristics.genotype. ))
bpr_blm <- filter(bpr_merge,bpr_merge$Characteristics.genotype.=="homozygous Blm knockout")  %>% select(-c(start,end,Characteristics.genotype. ))
bpr_blm$width=bpr_blm$CI.end-bpr_blm$CI.start
bpr_wt$width=bpr_wt$CI.end-bpr_wt$CI.start

bpr_blm=select(bpr_blm,c("seqnames"   ,"CI.start", "CI.end" , "genoT",   "Extract.Name", "width","filenames" ))
bpr_wt=select(bpr_wt,c("seqnames"  ,"CI.start","CI.end"  , "genoT",   "Extract.Name", "width" ,"filenames"))
colnames(bpr_blm) = c("chr"   ,  "start"     ,   "end"    , "genoT",   "Extract.Name", "width","filenames"    )
colnames(bpr_wt) = c("chr"   ,  "start"     ,   "end"  ,  "genoT",   "Extract.Name", "width"  ,"filenames"    )

bpr_blm <- select(bpr_blm,c(chr,start,end,width,genoT,Extract.Name,filenames))
bpr_wt <- select(bpr_wt,c(chr,start,end,width,genoT,Extract.Name,filenames))

bpr=rbind(bpr_blm,bpr_wt)
bpr$ploidy=2

write.table(bpr,"Precision Recall/Input/Niek_BPR_unfiltered.txt",quote = F,row.names = F,col.names = T,sep = "\t")

###############################################################################################################
##################### calculating precision recall for unfiltered SCEs from BPR 100reads ######################
###############################################################################################################

bpr=read.table("Precision Recall/Input/Niek_BPR_unfiltered.txt",header=T)

prec_recall=data.frame(matrix(0,nrow=3,ncol=5))
colnames(prec_recall)=c(0,10000,100000,500000,1e6)
n=1
sumOfTPs=0

for (file in levels(droplevels(bpr$Extract.Name))){
	message("Working on ",file," ... ",round((n/length(levels(droplevels(bpr$Extract.Name))))*100,digits = 1),"%")
	n=n+1
	tmp = GRanges(filter(bpr,bpr$Extract.Name==file))
	if(nrow(filter(bait,Extract.Name==file))!=0){
		tmp_bench= GRanges(filter(bait,Extract.Name==file))
		sumOfTPs=sumOfTPs+length(tmp_bench)
	} else{ 
		message("Not found",file)
		next
	}
	for (i in c(0,10000,100000,500000,1e6)){
		(TP <- length(countOverlaps(tmp_bench,tmp)[countOverlaps(tmp_bench,tmp,maxgap = i) !=0]))
		(FP=length(tmp)-TP)
		#(FP <-  length(tmp[-queryHits(findOverlaps(tmp_bench,tmp, type="any",maxgap = i)),]))
		(FN=length(tmp_bench)-TP)
		#(FN <- length(tmp_bench[-queryHits(findOverlaps(tmp_bench,tmp, type="any",maxgap = i)),]))
		prec_recall[1,as.character(i)]=prec_recall[1,as.character(i)]+TP
		prec_recall[2,as.character(i)]=prec_recall[2,as.character(i)]+FN
		prec_recall[3,as.character(i)]=prec_recall[3,as.character(i)]+FP
	}
}
pr = as.data.frame(t(prec_recall))
colnames(pr)=c("TP","FN","FP")
pr$precision=(pr$TP/(pr$TP+pr$FP))
pr$recall=(pr$TP/(pr$TP+pr$FN))
pr$gap=rownames(pr)
pr$gap<- as.factor(as.character(pr$gap))
pr$type ="unfiltered"
pr_unfilt=pr
ggplot(pr)+geom_point(aes(recall,precision,shape=gap,group=type,color=type),size=4)+
	geom_line(aes(recall,precision,group=type,color=type),size=1)+
	xlim(0,1)+ylim(0,1)+
	theme_bw()

#############################################################################################
##################### filtering Niek's SCEs from BPR using custom method #####################
##############################################################################################

bpr = read.table("Precision Recall/Input/Niek_BPR_unfiltered.txt",header=T)
fucci_breakpoints_good = bpr

# 1) Remove homozygous events
fucci_breakpoints_good_hetero = filter(fucci_breakpoints_good,(ploidy==1 & (genoT== "cc-ww"| genoT== "ww-cc" )) | ( ploidy==2 & (genoT!= "cc-ww"& genoT!= "ww-cc" )))
fucci_breakpoints_good_hetero$library<-as.factor(fucci_breakpoints_good_hetero$Extract.Name)

# 2) Filter out events that are too close to each other (2Mb)
fucci_breakpoints_good_hetero_SCE = data.frame()
n=1
fucci_breakpoints_good_hetero=select(fucci_breakpoints_good_hetero,c(chr,start,end,genoT,Extract.Name,ploidy,library))
colnames(fucci_breakpoints_good_hetero)=c("seqnames",  "start","end"  ,"genoT" ,"Extract.Name" ,"ploidy","library" )

for (level in levels(droplevels(fucci_breakpoints_good_hetero$library))){
	message(  round(     (      (n/length(levels(fucci_breakpoints_good_hetero$library)))   *100        )     ,digits = 2)   ,"% ... complete"   )
	tmp = filter(fucci_breakpoints_good_hetero, library==level)
	tmp$seqnames <- droplevels(tmp$seqnames)
	
	level2=levels(tmp$seqnames)[1]
	for (level2 in levels(tmp$seqnames)){
		
		tmp2 = filter(tmp, seqnames==level2)
		tmp3 <- GRanges(tmp2)
		
		overlaps = countOverlaps(tmp3, tmp3, type="any",maxgap = 2000000)
		overlap_rows = which(overlaps %in% c(1))
		
		tmp4 = tmp2[overlap_rows,]
		fucci_breakpoints_good_hetero_SCE <- rbind(fucci_breakpoints_good_hetero_SCE,tmp4)
		
	} n=n+1
}
fucci_breakpoints_good_hetero_SCE$library<-droplevels(fucci_breakpoints_good_hetero_SCE$library)

# 3) Filter out events near centromeres
centromeres <- read.table("/Users/zeidh/Desktop/THESIS/Anxilliary/centromeres2.txt",header=F,fill=T) #%>% select(-c(V4))
centromeres <-centromeres %>% dplyr::rename("seqnames"=V1,"start"=V2,"end"=V3)
centroGRange <- GRanges(centromeres)
summaryBreaks.df <- GRanges(fucci_breakpoints_good_hetero_SCE)
fucci_breakpoints_good_hetero_SCE_CM_bl <-as.data.frame(summaryBreaks.df[-queryHits(findOverlaps(summaryBreaks.df, centroGRange, type="any")),])
fucci_breakpoints_good_hetero_SCE_CM_bl$library<-droplevels(fucci_breakpoints_good_hetero_SCE_CM_bl$library)

# Write to enrichment folder for enrichment analysis
write.table(fucci_breakpoints_good_hetero_SCE_CM_bl, "/Users/zeidh/Desktop/THESIS/Precision Recall/Input/Niek_BPR_filtered.txt", col.names = T, row.names = F, quote = F, sep="\t")

#########################################################################################
##################### calculating precision recall of filtered data #####################
#########################################################################################


fucci_breakpoints_good_hetero_SCE_CM_bl=read.table("/Users/zeidh/Desktop/THESIS/Precision Recall/Input/Niek_BPR_filtered.txt", header = T)

prec_recall=data.frame(matrix(0,nrow=3,ncol=5))
colnames(prec_recall)=c(0,10000,100000,500000,1e6)
n=1
sumOfTPs=0

for (file in levels(droplevels(fucci_breakpoints_good_hetero_SCE_CM_bl$Extract.Name))){
	message("Working on ",file," ... ",round((n/length(levels(droplevels(fucci_breakpoints_good_hetero_SCE_CM_bl$Extract.Name))))*100,digits = 1),"%")
	n=n+1
	tmp = GRanges(filter(fucci_breakpoints_good_hetero_SCE_CM_bl,fucci_breakpoints_good_hetero_SCE_CM_bl$Extract.Name==file))
	if(nrow(filter(bait,Extract.Name==file))!=0){
		tmp_bench= GRanges(filter(bait,Extract.Name==file))
		sumOfTPs=sumOfTPs+length(tmp_bench)
	} else{ 
		message("Not found",file)
		next
	}
	for (i in c(0,10000,100000,500000,1e6)){
		(TP <- length(countOverlaps(tmp_bench,tmp)[countOverlaps(tmp_bench,tmp,maxgap = i) !=0]))
		(FP=length(tmp)-TP)
		#(FP <-  length(tmp[-queryHits(findOverlaps(tmp_bench,tmp, type="any",maxgap = i)),]))
		(FN=length(tmp_bench)-TP)
		#(FN <- length(tmp_bench[-queryHits(findOverlaps(tmp_bench,tmp, type="any",maxgap = i)),]))
		prec_recall[1,as.character(i)]=prec_recall[1,as.character(i)]+TP
		prec_recall[2,as.character(i)]=prec_recall[2,as.character(i)]+FN
		prec_recall[3,as.character(i)]=prec_recall[3,as.character(i)]+FP
		if (TP+FP != length(tmp)){
			print("something doesnt add up")
		}
		if (TP+FN != length(tmp_bench)){
			print("something doesnt add up x2")
		}
	}
	
}
sumOfTPs
pr = as.data.frame(t(prec_recall))
colnames(pr)=c("TP","FN","FP")
pr$precision=(pr$TP/(pr$TP+pr$FP))
pr$recall=(pr$TP/(pr$TP+pr$FN))
pr$gap=rownames(pr)
pr$gap<- as.factor(as.character(pr$gap))
pr$type ="filtered"
pr_filt=pr
ggplot(pr)+geom_point(aes(recall,precision,shape=gap,group=type,color=type),size=4)+
	geom_line(aes(recall,precision,group=type,color=type),size=1)+
	xlim(0,1)+ylim(0,1)+
	theme_bw()

pr_merge=rbind(pr_filt,pr_unfilt)
write.table(pr_merge,"Precision Recall/Input/PrecisionRecall_filt_unfilt_BPR_SCEs_niek.txt",quote = F,row.names = F,col.names = T,sep = "\t")


#########################################################################################
##################### calculating precision recall of filtered data #####################
#########################################################################################

pr_merge=read.table("Precision Recall/Input/PrecisionRecall_filt_unfilt_BPR_SCEs_niek.txt",header = T)

pr_merge$gap<-as.factor(pr_merge$gap)


ggplot(pr_merge)+geom_point(aes(recall,precision,shape=gap,group=type,color=type),size=4)+
	geom_line(aes(recall,precision,group=type,color=type),size=1)+
	xlim(0,1)+ylim(0,1)+
	labs(x="Recall",y="Precision")+
	theme_classic()+
	theme(text = element_text(size = 15),aspect.ratio = 1)+
	ggsave("Precision Recall/prec_recall.png")



