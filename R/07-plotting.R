#######################################################
####         Workflow Part 8: Re-plotting BPR       ###
#######################################################
library(tidyverse)
library(plyr)
library(ggpubr)



breakpointSummary <- read.table("SCEs/fucci_sces_complete.bed",header=T) #%>% select(-c(CI.start,CI.end,genoT))





########### PLOT #1: SCEs/CHR/LIB vs LENGTH ###############
###########################################################

breakpointSummary$gene <- as.factor(breakpointSummary$gene)
breakpointSummary$library <- as.factor(breakpointSummary$library)

#sces per chroomosome
suppressMessages(all <- as.data.frame(breakpointSummary %>% group_by(seqnames) %>% dplyr::summarize(ALL=n())))
suppressMessages(a <- as.data.frame(breakpointSummary %>% filter(gene=="BLM")  %>% group_by(seqnames) %>% dplyr::summarize(BLM=n())))
suppressMessages(b <- as.data.frame(breakpointSummary %>% filter(gene=="BLM/RECQL5")  %>% group_by(seqnames)%>% dplyr::summarize("BLM/RECQL5"= n())))
suppressMessages(c <- as.data.frame(breakpointSummary %>% filter(gene=="RECQL5")  %>% group_by(seqnames)%>% dplyr::summarize(RECQL5=n())))
suppressMessages(d <- as.data.frame(breakpointSummary %>% filter(gene=="WT")  %>% group_by(seqnames)%>% dplyr::summarize(WT=n())))
suppressMessages(e <- as.data.frame(breakpointSummary %>% filter(gene=="RECQL1")  %>% group_by(seqnames)%>% dplyr::summarize(RECQL1=n())))
suppressMessages(f <- as.data.frame(breakpointSummary %>% filter(gene=="WRN")  %>% group_by(seqnames)%>% dplyr::summarize(WRN=n())))
suppressMessages(g <- as.data.frame(breakpointSummary %>% filter(gene=="WRN/RECQL5")  %>% group_by(seqnames)%>% dplyr::summarize("WRN/RECQL5"=n())))
suppressMessages(h <- as.data.frame(breakpointSummary %>% filter(gene=="RTEL1")  %>% group_by(seqnames)%>% dplyr::summarize(RTEL1=n())))


byChr <- merge(a,b,by="seqnames")
byChr <- merge(byChr,c,by="seqnames",all=T)
byChr <- merge(byChr,d,by="seqnames",all=T)
byChr <- merge(byChr,e,by="seqnames",all=T)
byChr <- merge(byChr,f,by="seqnames",all=T)
byChr <- merge(byChr,g,by="seqnames",all=T)
byChr <- merge(byChr,h,by="seqnames",all=T)
byChr <- merge(byChr,all,by="seqnames",all=T)
byChr[is.na(byChr)] <- 0
write.table(byChr,"/Users/zeidh/Desktop/THESIS/Anxilliary/perChrom.txt",quote=F,row.names = F,col.names = T,sep="\t")

breakpointSummary$Library=breakpointSummary$library
numOfLibsPerGene <- data.frame(gene=character(),n=numeric())
a <- (as.data.frame(breakpointSummary %>% filter(gene=="BLM")  %>% droplevels()))# group_by(seqnames)%>% summarize(n()))
numOfLibsPerGene= add_row(numOfLibsPerGene, gene="BLM",n=as.numeric(length(levels(a$Library))))
b <- (as.data.frame(breakpointSummary %>% filter(gene=="BLM/RECQL5")  %>% droplevels()))# group_by(seqnames)%>% summarize(n()))
numOfLibsPerGene= add_row(numOfLibsPerGene, gene="BLM.RECQL5",n=as.numeric(length(levels(b$Library))))
c <- (as.data.frame(breakpointSummary %>% filter(gene=="RECQL1")  %>% droplevels()))# group_by(seqnames)%>% summarize(n()))
numOfLibsPerGene= add_row(numOfLibsPerGene, gene="RECQL1",n=as.numeric(length(levels(c$Library))))
d <- (as.data.frame(breakpointSummary %>% filter(gene=="WT")  %>% droplevels()))# group_by(seqnames)%>% summarize(n()))
numOfLibsPerGene= add_row(numOfLibsPerGene, gene="WT",n=as.numeric(length(levels(d$Library))))
e <- (as.data.frame(breakpointSummary %>% filter(gene=="RECQL5")  %>% droplevels()))# group_by(seqnames)%>% summarize(n()))
numOfLibsPerGene= add_row(numOfLibsPerGene, gene="RECQL5",n=as.numeric(length(levels(e$Library))))
f <- (as.data.frame(breakpointSummary %>% filter(gene=="RTEL1")  %>% droplevels()))# group_by(seqnames)%>% summarize(n()))
numOfLibsPerGene= add_row(numOfLibsPerGene, gene="RTEL1",n=as.numeric(length(levels(f$Library))))
g <- (as.data.frame(breakpointSummary %>% filter(gene=="WRN")  %>% droplevels()))# group_by(seqnames)%>% summarize(n()))
numOfLibsPerGene= add_row(numOfLibsPerGene, gene="WRN",n=as.numeric(length(levels(g$Library))))
h <- (as.data.frame(breakpointSummary %>% filter(gene=="WRN/RECQL5")  %>% droplevels()))# group_by(seqnames)%>% summarize(n()))
numOfLibsPerGene= add_row(numOfLibsPerGene, gene="WRN.RECQL5",n=as.numeric(length(levels(h$Library))))
numOfLibsPerGene= add_row(numOfLibsPerGene, gene="ALL",n=as.numeric(length(levels(breakpointSummary$Library))))
write.table(numOfLibsPerGene,"/Users/zeidh/Desktop/THESIS/Anxilliary/numOfLibsPerGene.txt",quote=F,row.names = F,col.names = F,sep="\t")


numOfLibsPerGene=read.table("/Users/zeidh/Desktop/THESIS/Anxilliary//numOfLibsPerGene.txt",header=F)
lengths <- read.table("/Users/zeidh/Desktop/THESIS/Anxilliary//hg38_chrLengths.txt",header=F)
colnames(lengths)=c("Chromosome","Start","Length")
lengths=select(lengths,c(Chromosome,Length))
lengths$Chromosome<-paste0("chr",lengths$Chromosome)
byChr<- read.table("/Users/zeidh/Desktop/THESIS/Anxilliary//perChrom.txt",header=T)

for (i in c("BLM","BLM.RECQL5", "RECQL1" ,"RECQL5","RTEL1","WRN","WRN.RECQL5","WT","ALL")){
	print(byChr[,i])
	for (j in 1:nrow(numOfLibsPerGene)){
		if (numOfLibsPerGene[j,1]==i){
			print(numOfLibsPerGene[j,2])
			byChr[,i] = byChr[,i]/ numOfLibsPerGene[j,2]
		}
	}
}


scePerChrPerGeneVsLength <- merge(lengths,byChr,by.x="Chromosome",by.y="seqnames")
tidy <- gather(scePerChrPerGeneVsLength, gene,sce_per_chr,BLM:ALL)
all <- filter(tidy,gene=="ALL")
tidy <- filter(tidy,gene!="ALL")
tidy$Length<- tidy$Length/1e+6

ggplot(tidy) +  geom_smooth(aes(Length,sce_per_chr,group=gene,color=gene),se=F,method="lm", na.rm=TRUE,size=2)+
	theme_classic(base_size = 20) +
	ylab("SCEs/CHR/LIB")+
	xlab("CHR LENGTH")+
	theme(legend.position = c(0.25, 0.8),legend.title = element_blank())+
	geom_point(aes(Length,sce_per_chr,group=gene),color="blue", na.rm=TRUE)

tidy$Length<-as.numeric(tidy$Length)
tidy$gene<-as.factor(tidy$gene)
tidy <- select(tidy,-Chromosome)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
					  conf.interval=.95, .drop=TRUE) {
	library(plyr)
	
	# New version of length which can handle NA's: if na.rm==T, don't count them
	length2 <- function (x, na.rm=FALSE) {
		if (na.rm) sum(!is.na(x))
		else       length(x)
	}
	
	# This does the summary. For each group's data frame, return a vector with
	# N, mean, and sd
	datac <- ddply(data, groupvars, .drop=.drop,
				   .fun = function(xx, col) {
				   	c(N    = length2(xx[[col]], na.rm=na.rm),
				   	  mean = mean   (xx[[col]], na.rm=na.rm),
				   	  sd   = sd     (xx[[col]], na.rm=na.rm)
				   	)
				   },
				   measurevar
	)
	
	# Rename the "mean" column    
	datac <- rename(datac, c("mean" = measurevar))
	
	datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
	
	# Confidence interval multiplier for standard error
	# Calculate t-statistic for confidence interval: 
	# e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
	ciMult <- qt(conf.interval/2 + .5, datac$N-1)
	datac$ci <- datac$se * ciMult
	
	return(datac)
}

tgc=summarySE(tidy, measurevar="sce_per_chr", groupvars=c("Length"))
tgc$Length=as.numeric(tgc$Length)

ggplot(tgc, aes(x=Length, y=sce_per_chr)) + 
	geom_errorbar(aes(ymin=sce_per_chr-se, ymax=sce_per_chr+se),width=10,size=1,color="grey") +
	geom_smooth(method="lm",se=F,size=3,color="black") +
	#geom_crossbar(aes(ymin=sce_per_chr-se, ymax=sce_per_chr+se))+
	theme_classic(base_size = 25) +
	ylab("SCEs / Chromosome")+
	xlab("Chromosome Length")+
	geom_point(size=4,color="#4075a8")+
	stat_cor(method = "pearson",p.accuracy = 0.001,size=7) +
	ggsave("Output/scePerChrPerGeneVsLength.png")






################## Plot #2 ##############################

suppressMessages(test<-as.data.frame(breakpointSummary %>%
									 	group_by(Library) %>%
									 	dplyr::summarize(n())))
ploidy = read.table("/Users/zeidh/Desktop/THESIS/Metrics/05-pe.se.metrics.good.quality.ploidy.txt",header=T)%>% select(c(file,ploidy))

test=merge(test,ploidy,by.x="Library",by.y="file")
test$norm=test$`n()`/test$ploidy
test$gene <- "gene"
test$library <- test$Library
suppressWarnings(test <- test %>% separate(Library, c("a","b","c","d","e","f"), "[_-]+"))

for (row in 1:nrow(test)){
	
	for (letter in c("a","b","c","d","e","f")){
		if (is.na(test[row,letter])!=T){
			if (test[row,letter]=="WT" | test[row,letter]=="wt"){
				test[row,"gene"]="WT"
			}
			else if (test[row,letter]=="blm" | test[row,letter]=="BLM"  | test[row,letter]=="blm1" ) {
				test[row,"gene"]="BLM"
			}
			else if (test[row,letter]=="wrn"  ) {
				if (test[row,"a"]=="recql5"){
					test[row,"gene"]="WRN/RECQL5"
				}
				else{
					test[row,"gene"]="WRN"
				}
			}
			
			else if (test[row,letter]=="RECQL5" | test[row,letter]=="recql5" | test[row,letter]=="RECQ5" | test[row,letter]=="recq5" ) {
				if (test[row,"gene"]=="BLM"){
					test[row,"gene"]="BLM/RECQL5"
				}
				else{
					test[row,"gene"]="RECQL5"
				}
			}
			else if (test[row,letter]=="blmrecq5"){
				test[row,"gene"]="BLM/RECQL5"
			}				
			else if (test[row,letter]=="RECQL1"| test[row,letter]=="recql1"){
				test[row,"gene"]="RECQL1"
			}
			else if ( test[row,letter]=="fucciwt"| test[row,letter]=="kbm7" | test[row,letter]=="wtfucci" ){
				test[row,"gene"]="WT"
			}
			else if (test[row,letter]=="RTEL"| test[row,letter]=="rtel"){
				test[row,"gene"]="RTEL1"
			}
		}
	}
	
}

test <- select(test,c("norm","gene"))
test<-dplyr::rename(test,c("sces"="n()"))
test=test[test$gene!="RTEL1",]
pval= data.frame(gene=c("BLM"  ,"BLM/RECQL5", "RECQL1", "RECQL5", "WRN" , "WRN/RECQL5", "WT"),
				 sig = c("***","***","ns","ns","***","*",""))

p=ggplot(test) + geom_jitter(aes(gene,norm),alpha=0.8)+ geom_boxplot(aes(gene,norm),color="#d92323",width=0.1,coef = 10) +
	theme_classic()+
	#geom_text(data=pval,aes(x=gene,y=(max(test$norm)+5),label=sig)) +
	theme(text=element_text(size=15),legend.position = "none")+
	labs(x="RecQ Knockout",y="SCEs/Library")+
	ggsave("Output/SCEs_per_libirary.png")


my_comparisons <- list( c("WT", "RECQL5"), c("WT", "BLM") ,c("WT", "RECQL1"),c("WT","RTEL1"),c("WT", "WRN"),c("WT", "WRN/RECQL5"), c("WT", "BLM/RECQL5"))


df_p_val <- test %>% rstatix::t_test( norm ~ gene, comparisons = list( c("BLM","WT"),c("BLM","BLM/RECQL5"),c("BLM/RECQL5","WT"),c("WT", "RECQL1"),c("WT","RECQL5"),c("WRN","WT"),c("WT","WRN/RECQL5"))) %>% rstatix::add_xy_position()
df_p_val$y.position = 50
p1 <- p + add_pvalue(df_p_val,label = "p.adj.signif",label.size=3,step.increase = 0.06,tip.length = 0)
p1
p1+ggsave("Output/SCEs_per_libirary.png")







breakpointSummary$width=breakpointSummary$end- breakpointSummary$start

suppressMessages(suppressWarnings(ggplot(breakpointSummary) + geom_smooth(aes(reads.per.mb.bpr, width),se=F) +
								  	#scale_y_log10() +
								  	theme_classic() +
								  	theme(text=element_text(size=20))+
								  	geom_hline(yintercept=19215, linetype="dashed", color = "red")+
								  	labs(x="Sequencing Effort (Reads/Mb)",y="SCE Breakpont Resoluton")+
								  	ggsave("Output/resolution_vs_depth.png")))


sce_summary=data.frame(gene=character(),SCE=numeric(),mean_resolution=numeric(),median_resolution=numeric())
b<- filter(breakpointSummary, gene=="BLM")
sce_summary<- add_row(sce_summary,gene="BLM",SCE=nrow(b),mean_resolution=mean(b$width),median_resolution=median(b$width))
r<- filter(breakpointSummary, gene=="RECQL5")
sce_summary<- add_row(sce_summary,gene="RECQL5",SCE=nrow(r),mean_resolution=mean(r$width),median_resolution=median(r$width))
br<- filter(breakpointSummary, gene=="BLM/RECQL5")
sce_summary<- add_row(sce_summary,gene="BLM/RECQL5",SCE=nrow(br),mean_resolution=mean(br$width),median_resolution=median(br$width))
w<- filter(breakpointSummary, gene=="WT")
sce_summary<- add_row(sce_summary,gene="WT",SCE=nrow(w),mean_resolution=mean(w$width),median_resolution=median(w$width))
write.table(sce_summary,"/Users/zeidh/Desktop/THESIS/Anxilliary//SCE_summary.txt",quote=F,row.names = F,col.names = T,sep="\t")




#plot here
ggplot(breakpointSummary) + stat_ecdf(aes(width,group=gene),size=0.3,color="#949494") +
	stat_ecdf(aes(width),size=2) +
	scale_x_log10() +
	theme_classic() +
	ylab("SCEs Mapped (%)") +
	xlab("Resolution") +
	annotation_logticks(sides = "b") +
	theme(text = element_text(size=20))+
	geom_density(aes(width),size=1,color="#565c7a")+
	guides(group="none")+
	geom_vline(xintercept=median(breakpointSummary$width), linetype="dashed", color = "red") +
	geom_text(aes(x=1000, label=paste0("Median\n",median(width)," bp"), y=0.8))+
	ggsave("Output/resolution.png")
				 
				 
				 