###################################################
				#Packages#
###################################################
library(ggplot2)
library(ggprism)
library(patchwork)
library(magrittr)
library(tidyverse)
library(dplyr)


###################################################
					#CNAs per cell#
################################################### 

bind=read.table("CNAs/Data/05-cnas_per_cell_datset_somatic.txt",header = T)
bind=bind[bind$gene!="RTEL1",]
bind=droplevels(bind)

################################################### 

cnas=read.table("CNAs/Data/02-somatic_cnas.txt",header=T)
cnas2=read.table("CNAs/Data/02-somatic_cnas_wt.txt",header=T)
cnas=rbind(cnas,cnas2)
labs= as.data.frame(cnas %>% group_by(type,gene) %>% dplyr::summarize(n()))
colnames(labs)=c("type", "gene","n")
labs$type=gsub(x=labs$type,pattern = "gain","Duplications")
labs$type=gsub(x=labs$type,pattern = "loss","Deletions")

################################################### 

p=ggplot(data = bind, aes(x = gene, y =cnas )) +
	geom_violin(aes(fill=type),scale = "count", trim = F) +
	geom_text(data=labs,aes(x=gene,y=300,label=paste0(n)),size=3)+
	theme_classic()+
	scale_y_log10()+
	geom_boxplot(width=0.1,fill="white")+
	theme(legend.position = "none",axis.title.x = element_blank(),text=element_text(size=15))+
	labs(x="Helicase KO",y="CNA / cell")+
	scale_fill_manual(values=c("#c96d6d", "#6d96c9"))+
	scale_x_discrete(breaks = c('BLM', 'BLM/RECQL5', 'RECQL1', 'RECQL5','WRN', 'WRN/RECQL5',"WT"),labels = c('BLM', 'BLM\nRECQL5', 'RECQL1', 'RECQL5','WRN', 'WRN\nRECQL5',"WT"))+
	facet_wrap(~type)
df_p_val <- bind %>% group_by(type) %>% rstatix::t_test( cnas ~ gene, comparisons = list( c("BLM","WT"),c("BLM/RECQL5","WT"),c("WT","RECQL5"),c("WRN","WT"),c("WT","WRN/RECQL5"))) %>% rstatix::add_xy_position()
df_p_val$y.position = c(1.3, 1.3 ,1.3 ,1.3 ,1.3, 0.7,0.7,0.7,0.7,0.7)# 1.3 1.3 1.3 1.3 1.3
p1 <- p + add_pvalue(df_p_val,label = "p.adj.signif",label.size=3,step.increase = 0.06,tip.length = 0)
p1
p1 + ggsave('Output/CNAP-per_cell.png')


###################################################
					# CNAs Width #
################################################### 


merge=read.table("CNAs/Data/02-somatic_cnas.txt",header = T)
merge2=read.table("CNAs/Data/02-somatic_cnas_wt.txt",header = T)
merge=rbind(merge,merge2)
merge$width=merge$width/1000000

###################################################

labs= as.data.frame(merge %>% group_by(type,gene) %>% dplyr::summarize(n()))
colnames(labs)=c("type", "gene","n")
merge$type=gsub(x = merge$type,pattern = "loss",replacement = "Deletions")
merge$type=gsub(x = merge$type,pattern = "gain",replacement = "Duplications")
labs$type=gsub(x = labs$type,pattern = "loss",replacement = "Deletions")
labs$type=gsub(x = labs$type,pattern = "gain",replacement = "Duplications")

###################################################

p=ggplot(data = merge, 
	   aes(x = gene, y =width )) +
	geom_violin(aes(fill=type),scale = "count", trim = F) +
	geom_text(data=labs,aes(x=gene,y=10000,label=paste0("n = ",n)),size=3)+
	scale_y_log10(labels = c("0","20","40", "80","160", "320","640","1280","2560"),breaks=c(0,20,40,80,160,320,640,1280,2560)) +
	#geom_dotplot(binaxis = "y",stackdir="center",dotsize=0.3)+
	geom_boxplot(width=0.2,fill="white")+
	theme_classic()+
	theme(legend.position = "none",axis.title.x = element_blank(),text=element_text(size=15))+
	labs(x="Helicase KO",y="CNV Length (Mb)")+
	scale_fill_manual(values=c("#c96d6d", "#6d96c9"))+
	scale_x_discrete(breaks = c('BLM', 'BLM/RECQL5', 'RECQL1' ,'RECQL5', 'WRN' ,'WRN/RECQL5', 'WT'), labels = c('BLM', 'BLM\nRECQL5', 'RECQL1' ,'RECQL5', 'WRN' ,'WRN\nRECQL5', 'WT'))+
	facet_wrap(~type)

df_p_val <- merge %>% group_by(type) %>% rstatix::t_test( width~ gene, comparisons = list(   c("BLM","WT"),c("BLM/RECQL5","WT"),c("WT","RECQL5"),c("WRN","WT"),c("WT","WRN/RECQL5")  )) %>% rstatix::add_xy_position()
df_p_val$y.position = c(2.8,2.8,2.8,2.8,2.8,1.7,1.7,1.7,1.7,1.7)# 1.3 1.3 1.3 1.3 1.3
p1
p1 <- p + add_pvalue(df_p_val,label = "p.adj.signif",label.size=3,step.increase = 0.06,tip.length = 0)
p1 + ggsave("Output/CNVlength.png")


