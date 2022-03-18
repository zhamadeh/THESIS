for (level in list.files("Intersecting_BED/")){

		basename = strsplit(x=level,split = "[.]")[[1]][1]
		null= read.table(paste0("EnrichmentPlots/Data/",basename,"_null.txt"),header=T,fill=T)
		true= read.table(paste0("EnrichmentPlots/Data/",basename,"_true.txt"),header=T,fill=T)
		pval= read.table(paste0("EnrichmentPlots/Data/",basename,"_pval.txt"),header=T,fill=T)
		
		null=null %>% mutate(color = ifelse(gene == "WT" |gene == "NA878", yes = "#89E088",no= "#88BCE0"))
		null=null[null$gene!="RTEL1",]
		true=true[true$gene!="RTEL1",]
		pval=pval[pval$gene!="RTEL1",]
		
		
		
		null$gene=factor(null$gene,levels = c("WT","WRN_RECQL5", "RECQL1", "WRN"  , "RECQL5"  ,"BLM"    ,   "BLM_RECQL5"   ,"NA878"     ) )
		ggplot(null)+geom_violin(aes(gene,enrich,fill=color),lwd=0.9,trim=F) +
			geom_boxplot(width=0.05,aes(gene,enrich),lwd=0.9) +
			theme_classic(base_size = 15) +
			geom_point(data=true,aes(gene,enrich),fill="red",colour="black",pch=21, size=5)+
			theme(legend.position = "none",aspect.ratio = 0.41)+
			labs(x="Cells",y="Enrichment",title = basename)+
			scale_x_discrete(breaks = c("WT","WRN_RECQL5", "RECQL1", "WRN"  , "RECQL5"  ,"BLM"    ,   "BLM_RECQL5"   ,"NA878"   ),labels = c("WT","WRN\nRECQL5", "RECQL1", "WRN"  , "RECQL5"  ,"BLM"    ,   "BLM\nRECQL5"   ,"NA878"   ))+
			geom_text(data=pval,aes(x=gene,y=(max(true$enrich,null$enrich)+0.2),label=sig))  +
			ggsave(paste0("Replotted/",basename,".pdf"))
	
}
# pdf()
# pdf()
# g
# dev.off()


for (level in list.files("Enrichment/ORCA/Plots/")){
	if (level  != "Data"){
		basename = strsplit(x = level,split = "[.]")[[1]][1]
		pval= read.table(paste0("Enrichment/ORCA/Plots/Data/",basename,"_pval.txt"),header=T,fill=T)
		#pval=pval[pval$gene!="RTEL1",]
		
		if ("***" %in% pval$sig){
			if (pval[pval$gene=="WT",]$sig=="ns"){
				interest = pval[pval$sig=="***",]
				if (nrow(interest)<2){
					message("Found *** with ",interest$gene," and ",basename)
				} else {
					message("Found *** with multiple KOs and",basename)
				}
			}
		}
	}
}
