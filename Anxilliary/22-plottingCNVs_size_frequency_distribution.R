
######################################################################
####   Plotting Karyogram   ####
######################################################################

################ Step 1: load genomic range data ####################

# Needs to be in dataframe format, can be points or ranges
dups=read.table("CNAs/Data/04-somatic_duplications.txt",header=T) %>% select(c(seqnames,start,end,width,gene))
dels=read.table("CNAs/Data/03-somatic_deletions.txt",header = T) %>% select(c(seqnames,start,end,width,gene))
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
labs <- aggregate(i~gene,merge,sum)
labs$width<-NA
labs
merge$width<-merge$width/1000000
ggplot(data = merge, 
       aes(x = gene, y =width )) +
  geom_violin(aes(),fill = "#a3a3a3",scale = "count", trim = F) +
  geom_text(data=labs,aes(x=gene,y=320,label=paste0("n = ",i)))+
  scale_y_log10(labels = c("0","20","40", "80","160", "320"),breaks=c(0,20,40,80,160,320)) +
  geom_dotplot(binaxis = "y",stackdir="center",dotsize=0.5)+
  geom_boxplot(width=0.1,fill="white")+
  theme(legend.position = "none")+
  labs(x="Helicase KO",y="CNV Length (Mb)")

merge %>% group_by(gene) %>% summarize(mean(width))

p=data.frame(ID=levels(read.table("Metrics/05-pe.se.metrics.good.quality.ploidy.txt",header=T)$library))
p$library <- p$ID
suppressWarnings(p <- p %>% separate(ID, c("a","b","c","d","e","f"), "[_-]+"))
p$gene <- "gene"
for (row in 1:nrow(p)){
  
  for (letter in c("a","b","c","d","e","f")){
    if (is.na(p[row,letter])!=T){
      if (p[row,letter]=="WT" | p[row,letter]=="wt"){
        p[row,"gene"]="WT"
      }
      else if (p[row,letter]=="blm" | p[row,letter]=="BLM" ) {
        p[row,"gene"]="BLM"
      }
      
      else if (p[row,letter]=="RECQL5" | p[row,letter]=="recql5" | p[row,letter]=="RECQ5" | p[row,letter]=="recq5" ) {
        if (p[row,"gene"]=="BLM"){
          p[row,"gene"]="BLM/RECQL5"
        }
        else{
          p[row,"gene"]="RECQL5"
        }
      }
    }
  }
  
}
p=select(p,-c(a,b,c,d,e,f))
p = as.data.frame(p %>% group_by(gene)%>% dplyr::summarize(n()))

m=merge %>% group_by(gene,type)%>% dplyr::summarize(n())
q=merge(m,p,by="gene")
q$norm = (q$`n().x`/q$`n().y`)*100



ggplot(q)+geom_col(aes(gene,norm,fill=type),color="black",width=0.5)+
  scale_y_continuous( expand = c(0, 0),limits = c(0,95))+
  geom_text(aes(x=gene,y=90,label=paste0("n = ",`n().y`)))+
  theme(legend.position = c(0.6,0.7),legend.title = element_blank(),axis.title.x = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  labs(y="CNV frequency (%)")+
  scale_fill_manual(values=c("#c96d6d", "#6d96c9"))#+
  #ggsave("OUTPUT/Plots/CNVs_frequency_dist.png")
q

q$i=1
labs <- aggregate(i~gene,q,sum)
labs$norm<-NA



                                  