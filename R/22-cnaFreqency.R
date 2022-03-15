######################################################################
####   Plotting Karyogram   ####
######################################################################

################ Step 1: load genomic range data ####################

# Needs to be in dataframe format, can be points or ranges
dups=read.table("CNAs/somatic_duplications.txt",header=T) %>% select(c(seqnames,start,end,width,gene))
dels=read.table("CNAs/somatic_deletions.txt",header = T) %>% select(c(seqnames,start,end,width,gene))
#wt=data.frame(seqnames=c("chr8","chr15"),start=c(1,60800000),end=c(139884560,89400000),gene=c("WT","WT"))
#wt$width= wt$end-wt$start
#dups=rbind(dups,wt)
dups$type="duplications"
dels$type="deletions"
merge= rbind(dups,dels)

ggplot(merge,aes(gene,width)) + geom_violin(trim=F,adjust=1)+scale_y_log10()+
	geom_point()
geom_text(data=summ,aes(x=gene,label=n))


require(scales)
merge$i=1
labs=as.data.frame(merge %>% group_by(gene,type)%>% dplyr::summarize(n()))
#labs <- aggregate(i~c(gene,type),merge,sum)
#labs$width<-NA
#labs
merge$width<-merge$width/1000000
ggplot(data = merge, 
	   aes(x = gene, y =width )) +
	geom_violin(aes(fill = type),scale = "count", trim = F) +
	geom_text(data=labs,aes(x=gene,y=320,label=paste0("n = ",`n()`)))+
	scale_y_log10(labels = c("0","20","40", "80","160", "320"),breaks=c(0,20,40,80,160,320)) +
	geom_dotplot(binaxis = "y",stackdir="center",dotsize=0.2)+
	geom_boxplot(width=0.1,fill="white")+
	theme_bw()+
	theme(legend.position = "none",text=element_text(size=15))+
	facet_wrap(~type)+
	labs(x="Helicase KO",y="CNV Length (Mb)")+
	ggsave("Output/cnaSizeDistribuytion.png")

 merge %>% group_by(gene) %>% summarize(mean(width))

p=data.frame(ID=levels(read.table("Metrics/05-pe.se.metrics.good.quality.ploidy.txt",header=T)$file))
p$library <- p$ID
suppressWarnings(p <- p %>% separate(library, c("a","b","c","d","e","f"), "[_-]+"))
p$gene <- "gene"
for (row in 1:nrow(p)){
	
	for (letter in c("a","b","c","d","e","f")){
		if (is.na(p[row,letter])!=T){
			if (p[row,letter]=="WT" | p[row,letter]=="wt"){
				p[row,"gene"]="WT"
			}
			else if (p[row,letter]=="blm" | p[row,letter]=="BLM"  | p[row,letter]=="blm1" ) {
				p[row,"gene"]="BLM"
			}
			else if (p[row,letter]=="wrn"  ) {
				if (p[row,"a"]=="recql5"){
					p[row,"gene"]="WRN/RECQL5"
				}
				else{
					p[row,"gene"]="WRN"
				}
			}
			
			else if (p[row,letter]=="RECQL5" | p[row,letter]=="recql5" | p[row,letter]=="RECQ5" | p[row,letter]=="recq5" ) {
				if (p[row,"gene"]=="BLM"){
					p[row,"gene"]="BLM/RECQL5"
				}
				else{
					p[row,"gene"]="RECQL5"
				}
			}
			else if (p[row,letter]=="blmrecq5"){
				p[row,"gene"]="BLM/RECQL5"
			}				
			else if (p[row,letter]=="RECQL1"| p[row,letter]=="recql1"){
				p[row,"gene"]="RECQL1"
			}
			else if ( p[row,letter]=="fucciwt"| p[row,letter]=="kbm7" | p[row,letter]=="wtfucci" ){
				p[row,"gene"]="WT"
			}
			else if (p[row,letter]=="RTEL"| p[row,letter]=="rtel"){
				p[row,"gene"]="RTEL1"
			}
		}
	}
	
}
p <- select(p,-c(a,b,c,d,e,f))
p$gene=as.factor(p$gene)

p = as.data.frame(p %>% group_by(gene)%>% dplyr::summarize(n()))

m=merge %>% group_by(gene,type)%>% dplyr::summarize(n())
q=merge(m,p,by="gene")
q$norm = (q$`n().x`/q$`n().y`)*100



ggplot(q)+geom_col(aes(gene,norm,fill=type),color="black",width=0.5)+
	scale_y_continuous( expand = c(0, 0),limits = c(0,95))+
	geom_text(aes(x=gene,y=90,label=paste0("n = ",`n().y`)))+
	theme_bw()+
	theme(legend.position = c(0.6,0.7),legend.title = element_blank(),axis.title.x = element_blank(),text=element_text(size=15))+
	labs(y="CNV frequency (%)")+
	scale_fill_manual(values=c("#c96d6d", "#6d96c9"))+
	ggsave("Output/CNVs_frequency_dist.png")



