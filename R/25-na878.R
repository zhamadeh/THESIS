na=read.table("../na878_sces.bed.txt")
na$V4=na$V3-na$V2
na=filter(na,V4<25000)
na$V1=gsub(x=na$V1,pattern = "chr",replacement = "")
na=sample_n(na,size = 1500)
write.table(na,"../NA878.bed",quote = F,row.names = F,col.names = F,sep = "\t")
