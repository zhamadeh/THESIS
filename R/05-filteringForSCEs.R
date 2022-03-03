metrics = read.table("/Users/zeidh/Desktop/THESIS/Metrics/05-pe.se.metrics.good.quality.ploidy.txt",header=T,fill=T) %>% filter(quality_manual=="G")

breakpoints = collectBreaksAllFiles(list.files("/Users/zeidh/Desktop/THESIS/BPR/data/",full.names = T))

fucci_breakpoints_good = merge(breakpoints,metrics,by.y="file",by.x="library")

######################## Filtering data ######################

# 1) Remove homozygous events

fucci_breakpoints_good_hetero = filter(fucci_breakpoints_good,(ploidy==1 & (genoT== "cc-ww"| genoT== "ww-cc" )) | ( ploidy==2 & (genoT!= "cc-ww"& genoT!= "ww-cc" )))
#levels(fucci_breakpoints_good_hetero$library) <- factor(fucci_breakpoints_good_hetero$library)
#length(levels(fucci_breakpoints_good_hetero$library))
fucci_breakpoints_good_hetero$library<-as.factor(fucci_breakpoints_good_hetero$library)
# 2) Filter out events that are too close to each other (2Mb)
fucci_breakpoints_good_hetero_SCE = data.frame()
n=1
level= levels(fucci_breakpoints_good_hetero$library)[1]
for (level in levels(fucci_breakpoints_good_hetero$library)){
	message(  round(     (      (n/length(levels(fucci_breakpoints_good_hetero$library)))   *100        )     ,digits = 2)   ,"% ... complete"   )
	tmp = filter(fucci_breakpoints_good_hetero, library==level)
	tmp$seqnames <- droplevels(tmp$seqnames)
	
	level2=levels(tmp$seqnames)[1]
	for (level2 in levels(tmp$seqnames)){
		
		tmp2 = filter(tmp, seqnames==level2)
	
		tmp3 <- GRanges(tmp2)
		
		overlaps = countOverlaps(tmp3, tmp3, type="any",maxgap = 2500000)
		overlap_rows = which(overlaps %in% c(1))
		
		tmp4 = tmp2[overlap_rows,]
		
		fucci_breakpoints_good_hetero_SCE <- rbind(fucci_breakpoints_good_hetero_SCE,tmp4)
		
	}
	n=n+1
}
fucci_breakpoints_good_hetero_SCE$library<-droplevels(fucci_breakpoints_good_hetero_SCE$library)
length(levels(fucci_breakpoints_good_hetero_SCE$library))

# 3) Filter out events near centromeres
centromeres <- read.table("/Users/zeidh/Desktop/THESIS/Anxilliary/centromeres2.txt",header=F,fill=T) #%>% select(-c(V4))
centromeres <-centromeres %>% dplyr::rename("seqnames"=V1,"start"=V2,"end"=V3)
centroGRange <- GRanges(centromeres)
summaryBreaks.df <- GRanges(fucci_breakpoints_good_hetero_SCE)
fucci_breakpoints_good_hetero_SCE_CM_bl <-as.data.frame(summaryBreaks.df[-queryHits(findOverlaps(summaryBreaks.df, centroGRange, type="any")),])
fucci_breakpoints_good_hetero_SCE_CM_bl$library<-droplevels(fucci_breakpoints_good_hetero_SCE_CM_bl$library)

# 4) Filter out events that are too close to each other (2Mb)
fucci_breakpoints_good_hetero_SCE_CM_bl_filtered_10kb <- filter(fucci_breakpoints_good_hetero_SCE_CM_bl, width < 10000) # 3k 
fucci_breakpoints_good_hetero_SCE_CM_bl_filtered_10kb$library<-droplevels(fucci_breakpoints_good_hetero_SCE_CM_bl_filtered_10kb$library)


######################## Exporting Data ######################

# Write to enrichment folder for enrichment analysis
write.table(fucci_breakpoints_good_hetero_SCE_CM_bl, "/Users/zeidh/Desktop/THESIS/SCEs/fucci_sces_complete.bed", col.names = T, row.names = F, quote = F, sep="\t")


######################## Plotting Data ######################
# Take putative SCEs and reprint breaksPlot with just them


cleanRData <- function(data=fucci_breakpoints_good_hetero_SCE_CM_bl){
	n=1
	for (file in levels(data$library)){
		
		message(  round(     (      (n/length(levels(data$library)))   *100        )     ,digits = 2)   ,"% ... complete"   )
		filename = paste0("/Users/zeidh/Desktop/THESIS/BPR/data/",file,".RData")
		tmp = get(load(filename))
		
		seqinfo <- tmp$breaks@seqinfo
		seqnameLevels <- levels(tmp$breaks@seqnames)
		
		tmp$breaks <-GRanges(filter(data,library==tmp$ID) %>% select(c(seqnames,start,end,width,strand,genoT,deltaW)))
		tmp$breaks@seqinfo <- seqinfo
		levels(tmp$breaks@seqnames) <- seqnameLevels
		
		tmp$confint = tmp$confint[queryHits(findOverlaps(tmp$confint, tmp$breaks, type="any")),]
		#$confint<-GRanges(filter(breaks,library==tmp$ID) %>% select(c(seqnames,CI.start,CI.end,width,strand,genoT,deltaW)))
		tmp$confint@seqinfo <- seqinfo
		levels(tmp$confint@seqnames) <- seqnameLevels
		
		save(tmp, file=paste0("/Users/zeidh/Desktop/THESIS/BPR/data_bl/",file,"filtered.RData"))
		
		n=n+1
	}
	
}

cleanRData(fucci_breakpoints_good_hetero)
breakpointR::plotBreakpoints(files2plot = list.files("../fucci/data_filtered/",full.names = T),file = "breaksPlot_na12878_good7_hetero.pdf")

cleanRData(fucci_breakpoints_good_hetero_SCE)
breakpointR::plotBreakpoints(files2plot = list.files("../fucci/data_filtered/",full.names = T),file = "breaksPlot_na12878_good7_hetero_sce.pdf")

cleanRData(fucci_breakpoints_good_hetero_SCE_CM_bl)
breakpointR::plotBreakpoints(files2plot = list.files("/Users/zeidh/Desktop/THESIS/BPR/data_bl/",full.names = T),file = "/Users/zeidh/Desktop/THESIS/breaksPlot_na12878_good7_hetero_sce_cm.pdf")

cleanRData(fucci_breakpoints_good_hetero_SCE_CM_bl_filtered_10kb)
breakpointR::plotBreakpoints(files2plot = list.files("../fucci/data_filtered/",full.names = T),file = "breaksPlot_na12878_good7_SCE_CM_10KB.pdf")