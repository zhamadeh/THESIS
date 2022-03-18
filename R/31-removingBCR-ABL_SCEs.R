############################################################################
############################	 Upload SCEs	##########################
############################################################################
sces=read.table("SCEs/fucci_sces_complete.bed",header=T)

############################################################################
############################   separate chr 9/22	########################
############################################################################
nine=sces[sces$seqnames=="chr9",]
nine=GRanges(nine)
twenty2=sces[sces$seqnames=="chr22",]
twenty2=GRanges(twenty2)

############################################################################
###########################   BCR/ABL coordinates ##########################
############################################################################

bcr=data.frame(seqnames=c("chr9","chr22"),start=c(130713015,23180364),end=c(130885683,23315950))
bcr=GRanges(bcr)

############################################################################
############################	 Check lengths	##########################
############################################################################

length(bcr)
length(nine)
length(twenty2)
length(subsetByOverlaps(nine,bcr))
length(subsetByOverlaps(nine,bcr,invert = T))
length(subsetByOverlaps(twenty2,bcr))
length(subsetByOverlaps(twenty2,bcr,invert = T))

############################################################################
###################### SCEs overlapping with gene fusion ####################
############################################################################
removeNine=subsetByOverlaps(nine,bcr)
removeTwenty2=subsetByOverlaps(twenty2,bcr)

sces=GRanges(sces)
sces=subsetByOverlaps(sces,removeNine,invert = T)
sces=subsetByOverlaps(sces,removeTwenty2,invert = T)

############################################################################
write.table(sces,"SCEs/sces_complete_minisPhili.txt",quote = F,row.names = F,col.names = T,sep = "\t")
############################################################################

############################################################################
###################### Central breakpoint for chr22 ####################
############################################################################

philli=GRanges(seqnames = "chr22",ranges = c(23289624:23291974))
potentialSCEs_22=subsetByOverlaps(removeTwenty2,philli,invert = T) # Look at non-intersecting with central breakppint
tt=dplyr::select(as.data.frame(potentialSCEs_22),c(seqnames,start,end,library))
write.table(tt,"../twentyTwoPossibleSCEs_2.bed",quote = F,row.names = F,col.names = F,sep = "\t")
# Export file names for concatenation
readsFiles=levels(droplevels(as.factor(tt$library)))
write.table(paste0(readsFiles,"_reads.bed.gz"),"../twentyTwoPossibleSCEs_files.txt",quote = F,row.names = F,col.names = F,sep = "\t")

# xargs -a twentyTwoPossibleSCEs_files.txt cp -t TwentyTwo/
# cat TwentyTwo/* > TwentyTwo/cat_twentytwo_reads.bed.gz
# gzip -d cat_twentytwo_reads.bed.gz
# grep "chr22\|track" cat_twentytwo_reads.bed >  cat_22_reads.bed


write.table(removeNine,"../removeNine.bed",quote = F,row.names = F,col.names = F,sep = "\t")

philli=GRanges(seqnames = "chr9",ranges = c(23289624:23291974))
potentialSCEs_22=subsetByOverlaps(removeTwenty2,philli,invert = T) # Look at non-intersecting with central breakppint
tt=dplyr::select(as.data.frame(potentialSCEs_22),c(seqnames,start,end,library))
write.table(tt,"../twentyTwoPossibleSCEs_2.bed",quote = F,row.names = F,col.names = F,sep = "\t")
# Export file names for concatenation
readsFiles=levels(droplevels(as.factor(tt$library)))
write.table(paste0(readsFiles,"_reads.bed.gz"),"../twentyTwoPossibleSCEs_files.txt",quote = F,row.names = F,col.names = F,sep = "\t")

# xargs -a twentyTwoPossibleSCEs_files.txt cp -t TwentyTwo/
# cat TwentyTwo/* > TwentyTwo/cat_twentytwo_reads.bed.gz
# gzip -d cat_twentytwo_reads.bed.gz
# grep "chr22\|track" cat_twentytwo_reads.bed >  cat_22_reads.bed








tt=dplyr::select(as.data.frame(removeTwenty2),c(seqnames,start,end,library))
tt=tt[order(tt$start),]
start=tt[1:5,]
end=tt[111:115,]
write.table(start,"../twentyTwo_left.bed",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(end,"../twentyTwo_right.bed",quote = F,row.names = F,col.names = F,sep = "\t")
# Export file names for concatenation
leftFiles=levels(droplevels(as.factor(start$library)))
write.table(paste0(leftFiles,"_reads.bed.gz"),"../twentyTwo_leftFiles.txt",quote = F,row.names = F,col.names = F,sep = "\t")
rightFiles=levels(droplevels(as.factor(end$library)))
write.table(paste0(rightFiles,"_reads.bed.gz"),"../twentyTwo_rightFiles.txt",quote = F,row.names = F,col.names = F,sep = "\t")



