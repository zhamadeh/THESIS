
m=read.table("Enrichment/Genes/essentialGenes_KBM.bed",header=F)
colnames(m)=c( "seqnames" ,"start" ,"end" ,"V4","V5")
m=m[m$seqnames %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),]
m$seqnames=paste0("chr",m$seqnames)
ESS=GRanges(m)

m=read.table("Enrichment/Genes/minusEssentialGenes_KBM.bed",header=F)
colnames(m)=c( "seqnames" ,"start" ,"end" ,"V4","V5")
m=m[m$seqnames %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),]
m$seqnames=paste0("chr",m$seqnames)
ESSminus=GRanges(m)

m=read.table("Enrichment/Genes/plusEssentialGenes_KBM.bed",header=F)
colnames(m)=c( "seqnames" ,"start" ,"end" ,"V4","V5")
m=m[m$seqnames %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),]
m$seqnames=paste0("chr",m$seqnames)
ESSplus=GRanges(m)


m=read.table("Enrichment/TRC/CD_HO.bed",header=F)
colnames(m)=c( "seqnames" ,"start" ,"end" ,"V4")
m=m[m$seqnames %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),]
m$seqnames=paste0("chr",m$seqnames)
CD_HO=GRanges(m)

m=read.table("Enrichment/TRC/CD_CD.bed",header=F)
colnames(m)=c( "seqnames" ,"start" ,"end" ,"V4")
m=m[m$seqnames %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),]
m$seqnames=paste0("chr",m$seqnames)
CD_CD=GRanges(m)
m=read.table("Enrichment/TRC/HO_HO.bed",header=F)
colnames(m)=c( "seqnames" ,"start" ,"end" ,"V4")
m=m[m$seqnames %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),]
m$seqnames=paste0("chr",m$seqnames)
HO_HO=GRanges(m)
m=read.table("Enrichment/TRC/HO_CD.bed",header=F)
colnames(m)=c( "seqnames" ,"start" ,"end" ,"V4")
m=m[m$seqnames %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),]
m$seqnames=paste0("chr",m$seqnames)
HO_CD=GRanges(m)







levels=list.files("Enrichment/TRC/",full.names = T)
geneSet=c("ESS","ESSplus","ESSminus")[1]
for (geneSet in c("ESS","ESSplus","ESSminus")){
	B=get(geneSet)
	
	null=data.frame()
	true=data.frame()
	pval=data.frame()
	for (level in levels){
		print(level)
		
		m=read.table(level,header=F)
		colnames(m)=c( "seqnames" ,"start" ,"end" ,"V4")
		m=m[m$seqnames %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),]
		m$seqnames=paste0("chr",m$seqnames)
		A=GRanges(m)
		
		pt <- overlapPermTest(A=A, B=B, ntimes=1000,mc.set.seed=FALSE,force.parallel=TRUE,mask=mask)
		tmp1=data.frame(Overlaps=pt$numOverlaps$permuted,gene=level)
		tmp1$enrich = tmp1$Overlaps/mean(tmp1$Overlaps)
		tmp2=data.frame(Overlaps=pt$numOverlaps$observed,gene=level)
		tmp2$enrich=tmp2$Overlaps/mean(tmp1$Overlaps)
		null=rbind(tmp1,null)
		true=rbind(tmp2,true)
		p = data.frame(gene=level,pval=pt$numOverlaps$pval)
		pval=rbind(pval,p)
	}
	pval=pval %>% mutate(sig = case_when(pval <= 0.001 ~ '***',pval < 0.01 ~ '**',pval < 0.05 ~ '*',pval < 1 ~ 'ns'))
	
	########################################################################################
	################################# Plotting #############################################
	########################################################################################
	
	null=null %>% mutate(color = ifelse(gene == "WT", "#89E088", "#88BCE0"))

	ggplot(null)+geom_violin(aes(gene,enrich,fill=color),lwd=0.9,trim=F) +
		geom_boxplot(width=0.05,aes(gene,enrich),lwd=0.9) +
		theme_classic(base_size = 15) +
		geom_point(data=true,aes(gene,enrich),fill="red",colour="black",pch=21, size=5)+
		theme(legend.position = "none",aspect.ratio = 0.41)  +
		labs(x="Cells",y="Enrichment",title=geneSet)+
		geom_text(data=pval,aes(x=gene,y=(max(true$enrich,null$enrich)+0.2),label=sig)) 
}
