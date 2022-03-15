sces= read.table("../fucci_sces.bed",header=T)%>%select(c(seqnames,start,end,gene,width))
sces$gene=gsub(x = sces$gene,pattern = "/",replacement = "_")
sces$gene=as.factor(sces$gene)
sces$seqnames=gsub(x = sces$seqnames,pattern = "chr",replacement = "")

for (level in levels(sces$gene)){
	print(level)
	tmp = filter(sces,gene==level)
	
	tmp = filter(tmp,width<25000)
	write.table(select(tmp,c(seqnames,start,end,gene)),paste0("../Enrichment/SCEs/",level,".bed"),quote = F,row.names = F,col.names = F,sep = "\t")
}

nrow(filter(sces,width<50000))

