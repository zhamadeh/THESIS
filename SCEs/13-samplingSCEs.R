


sces= read.table("SCEs/fucci_sces_complete.bed",header=T)%>%select(c(seqnames,start,end,gene,width,library))
sces$gene=gsub(x = sces$gene,pattern = "/",replacement = "_")
sces$gene=as.factor(sces$gene)

#sces = filter(sces,width<10000)

sample = sample(sces$library,size=10)

sampleChr = sample(levels(sces$seqnames),size=10)


breakpointR::plotBreakpointsPerChr(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data_bl/",sample,"filtered.RData"),plotspath = "/Users/zeidh/Desktop/THESIS/SCEs/Samples/Sample/",chromosomes = sampleChr)
breakpointR::plotBreakpointsPerChr(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data/",sample,".RData"),plotspath = "/Users/zeidh/Desktop/THESIS/SCEs/Samples/Sample_bl/",chromosomes = sampleChr)

breakpointR::plotBreakpointsPerChr(files2plot = paste0("/Users/zeidh/Desktop/THESIS/BPR/data_bl/",levels(droplevels(WT$library)),"filtered.RData"),
								   plotspath = "/Users/zeidh/Desktop/THESIS/SCEs/Sample/",chromosomes = levels(droplevels(WT$seqnames)))

breakpointR::plotBreakpoints(files2plot = paste0("BPR/data_bl/",sample(levels(droplevels(WT$library)),size=15 ),"filtered.RData"),file = "wt_sample.pdf")
