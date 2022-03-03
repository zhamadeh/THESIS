###################################################################################
##################     Workflow Part 1: Collect initial metrics   #################
###################################################################################

bammetrics=read.table("Desktop/THESIS/Metrics/thesis.bam.metrics.txt",header=T)
bammetrics_se=read.table("Desktop/THESIS/Metrics/thesis.bam.metrics_se.txt",header=T)
bammetrics$end = "pe"
bammetrics_se$end = "se"
`%notin%` <- Negate(`%in%`)
bammetrics=bammetrics[bammetrics$file %notin% bammetrics_se$file,]
bammetrics=rbind(bammetrics,bammetrics_se)


bprmetrics=read.table("Desktop/THESIS/Metrics/thesis.bpr.lib.metrics.txt",header=T)
bprmetrics_se=read.table("Desktop/THESIS/Metrics/thesis.bpr.metrics.se.txt",header=T)
bprmetrics$end = "pe"
bprmetrics_se$end = "se"
colnames(bprmetrics)=c("library","background" , "reads.per.mb.bpr" ,"coverage.bpr" ,"end" )
colnames(bprmetrics_se)=c("library","background" , "reads.per.mb.bpr" ,"coverage.bpr" ,"end" )
bprmetrics=bprmetrics[bprmetrics$library %notin% bprmetrics_se$library,]
bprmetrics=rbind(bprmetrics,bprmetrics_se)



merge = merge(bammetrics,bprmetrics,by.x="file",by.y="library")


merge$library <- merge$file
suppressWarnings(merge <- merge %>% separate(library, c("a","b","c","d","e","f"), "[_-]+"))
merge$gene <- "gene"
for (row in 1:nrow(merge)){
	
	for (letter in c("a","b","c","d","e","f")){
		if (is.na(merge[row,letter])!=T){
			if (merge[row,letter]=="WT" | merge[row,letter]=="wt"){
				merge[row,"gene"]="WT"
			}
			else if (merge[row,letter]=="blm" | merge[row,letter]=="BLM"  | merge[row,letter]=="blm1" ) {
				merge[row,"gene"]="BLM"
			}
			else if (merge[row,letter]=="wrn"  ) {
				if (merge[row,"a"]=="recql5"){
					merge[row,"gene"]="WRN/RECQL5"
				}
				else{
					merge[row,"gene"]="WRN"
				}
			}
			
			else if (merge[row,letter]=="RECQL5" | merge[row,letter]=="recql5" | merge[row,letter]=="RECQ5" | merge[row,letter]=="recq5" ) {
				if (merge[row,"gene"]=="BLM"){
					merge[row,"gene"]="BLM/RECQL5"
				}
				else{
					merge[row,"gene"]="RECQL5"
				}
			}
			else if (merge[row,letter]=="blmrecq5"){
				merge[row,"gene"]="BLM/RECQL5"
			}				
			else if (merge[row,letter]=="RECQL1"| merge[row,letter]=="recql1"){
				merge[row,"gene"]="RECQL1"
			}
			else if ( merge[row,letter]=="fucciwt"| merge[row,letter]=="kbm7" | merge[row,letter]=="wtfucci" ){
				merge[row,"gene"]="WT"
			}
			else if (merge[row,letter]=="RTEL"| merge[row,letter]=="rtel"){
				merge[row,"gene"]="RTEL1"
			}
		}
	}
	
}
merge <- select(merge,-c(a,b,c,d,e,f))
merge=merge[merge$gene!="gene",]
lessThan30=filter(merge,merge$reads.per.mb.bpr<=30)
merge=filter(merge,merge$reads.per.mb.bpr>30)
merge%>%group_by(gene)%>% summarize(n())

write.table(filter(merge,end.x=="se")$file,"Desktop/THESIS/se.list.of.libaries.over.30.rpm.txt",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(filter(merge,end.x=="pe")$file,"Desktop/THESIS/pe.list.of.libaries.over.30.rpm.txt",row.names = F,col.names = F,quote = F,sep = "\t")

write.table(filter(merge,end.x=="se"),"Desktop/THESIS/se.metrics.30.rpm.txt",row.names = F,col.names = T,quote = F,sep = "\t")
write.table(filter(merge,end.x=="pe"),"Desktop/THESIS/pe.metrics.30.rpm.txt",row.names = F,col.names = T,quote = F,sep = "\t")

write.table(filter(lessThan30,end.x=="se"),"Desktop/THESIS/se.metrics.0.rpm.txt",row.names = F,col.names = T,quote = F,sep = "\t")
write.table(filter(lessThan30,end.x=="pe"),"Desktop/THESIS/pe.metrics.0.rpm.txt",row.names = F,col.names = T,quote = F,sep = "\t")

write.table(filter(  (filter(pair,end=="se"))  , quality_manual == "G"  )$file,"Desktop/THESIS/se.good.libraries.txt",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(filter(  (filter(pair,end=="pe"))  , quality_manual == "G"  )$file,"Desktop/THESIS/pe.good.libraries.txt",row.names = F,col.names = F,quote = F,sep = "\t")





