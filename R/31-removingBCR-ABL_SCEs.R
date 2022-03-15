sces=read.table("SCEs/fucci_sces_complete.bed",header=T)
nine=sces[sces$seqnames=="chr9",]
nine=GRanges(nine)
twenty2=sces[sces$seqnames=="chr22",]
twenty2=GRanges(twenty2)

bcr=data.frame(seqnames=c("chr9","chr22"),start=c(130713015,23180364),end=c(130885683,23315950))
bcr=GRanges(bcr)

sces=GRanges(sces)

length(bcr)
length(nine)
length(twenty2)
length(subsetByOverlaps(nine,bcr))
length(subsetByOverlaps(nine,bcr,invert = T))
length(subsetByOverlaps(twenty2,bcr))
length(subsetByOverlaps(twenty2,bcr,invert = T))
removeNine=subsetByOverlaps(nine,bcr)
removeTwenty2=subsetByOverlaps(twenty2,bcr)
length(subsetByOverlaps(nine,bcr,invert = T))
length(subsetByOverlaps(twenty2,bcr,invert = T))

length(sces)
sces=subsetByOverlaps(sces,removeNine,invert = T)
sces=subsetByOverlaps(sces,removeTwenty2,invert = T)
length(sces)
write.table(sces,"SCEs/sces_complete_minisPhili.txt",quote = F,row.names = F,col.names = T,sep = "\t")

