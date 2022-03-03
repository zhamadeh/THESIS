########################################################################################
################## Plotting Niek's BS data with Victor's enrichment ####################
########################################################################################

########################################################################################
################################### Assembling metadata ##############################
########################################################################################

metadata <- read.table("Anxilliary/Metadata/E-MTAB-5976.sdrf.txt",sep="\t",header=T)
metadata <-  select(metadata,c(Comment.ENA_RUN.,Characteristics.cell.line., Extract.Name,Factor.Value.protocol., Characteristics.genotype., Characteristics.organism.))

#  Put all bed files in one file
bed <- data.frame()
for (file in list.files("Anxilliary/Benchmark_sce_files/",full.names = T)){
	tmp <- read.table(file,sep="\t")
	if (ncol(tmp)==4){
		tmp$V5=tmp$V3 - tmp$V2
		bed <- rbind(tmp,bed)
	} else {bed <- rbind(tmp,bed)}
	#message("columns: ",ncol(tmp)," row namese  ", nrow(tmp) )
}
#merge with metadata
merge <- merge(bed,metadata,by.x="V4",by.y="Extract.Name")

#sort into BLM and WT SCEs
wt <- filter(merge,merge$Characteristics.genotype.=="wild type genotype") %>% select(-c(Characteristics.organism., Factor.Value.protocol.,Comment.ENA_RUN.)) %>% select(c(V1,V2,V3,V5,V4,Characteristics.cell.line.,Characteristics.genotype.))
blm <- filter(merge,merge$Characteristics.genotype.=="homozygous Blm knockout")%>% select(-c(Characteristics.organism., Factor.Value.protocol.,Comment.ENA_RUN.))%>% select(c(V1,V2,V3,V5,V4,Characteristics.cell.line.,Characteristics.genotype.))
colnames(blm) = c("chr"   ,  "start"     ,   "end"  , "width"    ,    "Extract.Name" ,"Characteristics.cell.line.","Characteristics.genotype.")
blm <- blm %>% select(c(chr,start,end,width,Extract.Name,Characteristics.cell.line.,Characteristics.genotype.))
colnames(wt) = c("chr"   ,  "start"     ,   "end"  , "width"    ,    "Extract.Name","Characteristics.cell.line.","Characteristics.genotype." )
wt <- wt %>% select(c(chr,start,end,width,Extract.Name,Characteristics.cell.line.,Characteristics.genotype.))

niek=rbind(blm,wt)
niek= filter(niek,width<10000)

niek$Characteristics.cell.line.=gsub(x = niek$Characteristics.cell.line.,pattern = "GM02085",replacement = "BS1")
niek$Characteristics.cell.line.=gsub(x = niek$Characteristics.cell.line.,pattern = "GM03402",replacement = "BS2")
niek$Characteristics.cell.line.=gsub(x = niek$Characteristics.cell.line.,pattern = "GM16375",replacement = "BS3")
niek$Characteristics.cell.line.=gsub(x = niek$Characteristics.cell.line.,pattern = "GM17361",replacement = "BS4")

niek$Characteristics.cell.line.=gsub(x = niek$Characteristics.cell.line.,pattern = "GM07492",replacement = "WT1")
niek$Characteristics.cell.line.=gsub(x = niek$Characteristics.cell.line.,pattern = "GM07545",replacement = "WT2")
niek$Characteristics.cell.line.=gsub(x = niek$Characteristics.cell.line.,pattern = "GM12891",replacement = "WT3")
niek$Characteristics.cell.line.=gsub(x = niek$Characteristics.cell.line.,pattern = "GM12892",replacement = "WT4")


for (level in levels(droplevels(as.factor(niek$Characteristics.cell.line.)))){
	tmp=filter(niek,Characteristics.cell.line.==level)
	
	tmp$chr=gsub(x = tmp$chr,pattern = "chr",replacement = "")
	write.table(tmp,paste0("Enrichment/GenomePer/Niek/Niek_SCE_10Kb//",level,".bed"),quote = F,row.names = F,col.names = F,sep = "\t")
}

####################################################################################
############################## G4s HG37  ##########################################
####################################################################################

# Need to run Victor's enichment perl script and store output in "Enrichment/GenomePer/Niek/___/" 

####################################################################################
############################## G4s HG37  ##########################################
####################################################################################

null_full=data.frame()
true_full=data.frame()
pval=data.frame()

for (file in list.files("Enrichment/GenomePer/Niek/Niek_G4_hg37/",full.names = T)){
	tmp=read.table(file,header=T)
	true = data.frame(Shift=tmp[1,1],Overlaps=tmp[1,2],gene=strsplit(x=basename(file),split = "[.]")[[1]][1])
	null=tmp[2:nrow(tmp),]
	null$gene=strsplit(x=basename(file),split = "[.]")[[1]][1]
	
	true$enrich = true$Overlaps/mean(null$Overlaps)
	null$enrich = null$Overlaps/mean(null$Overlaps)
	
	null_full=rbind(null_full,null)
	true_full=rbind(true_full,true)
}
	
pval=data.frame(File=c("BS1","BS2","BS3","BS4","WT1","WT2","WT3","WT4"), P_value=c(0.012,0.001,0.235,0.001,0.082,0.032,0.257,0.347))

pval=pval %>%mutate(sig = case_when(P_value <= 0.001 ~ '***',P_value < 0.01 ~ '**', P_value < 0.05 ~ '*',P_value < 1 ~ 'ns'))
	
null_full=null_full %>% mutate(color = ifelse(gene %in% c("WT1","WT2","WT3","WT4") , "#89E088", "#88BCE0"))
null_full$gene <- factor(null_full$gene, levels = c( "WT1" ,"WT2" ,"WT3", "WT4","BS1" ,"BS2" ,"BS3", "BS4"))

ggplot(null_full)+geom_violin(aes(gene,enrich,fill=color),lwd=0.9,trim=F) +
	geom_boxplot(width=0.05,aes(gene,enrich),lwd=0.9) +
	#scale_fill_manual(values=c("#53b0db","#53b0db","#53b0db","#6ceb70"))+
	theme_classic(base_size = 15) +
	geom_point(data=true_full,aes(gene,enrich),fill="red",colour="black",pch=21, size=5)+
	theme(legend.position = "none")  +
	geom_text(data=pval,aes(x=File,y=(max(true_full$enrich,null_full$enrich)+0.2),label=sig))  +
	labs(x="Cells",y="Enrichment",title="G4s GRCh37") +
	ggsave("Enrichment/GenomePer/Niek/Niek_Plots/G4s_GRCh37.png")
	
	
################################################################################
############################## Genes HG37 ######################################
################################################################################
	
null_full=data.frame()
true_full=data.frame()
pval=data.frame()

for (file in list.files("Enrichment/GenomePer/Niek/Niek_Genes_hg37///",full.names = T)){
	tmp=read.table(file,header=T)
	true = data.frame(Shift=tmp[1,1],Overlaps=tmp[1,2],gene=strsplit(x=basename(file),split = "[.]")[[1]][1])
	null=tmp[2:nrow(tmp),]
	null$gene=strsplit(x=basename(file),split = "[.]")[[1]][1]
	
	true$enrich = true$Overlaps/mean(null$Overlaps)
	null$enrich = null$Overlaps/mean(null$Overlaps)
	
	null_full=rbind(null_full,null)
	true_full=rbind(true_full,true)
}

pval=data.frame(File=c("BS1","BS2","BS3","BS4","WT1","WT2","WT3","WT4"), P_value=c(0.001,0.001,0.151,0.001,0.083,0.034,0.306,0.407))

pval=pval %>%mutate(sig = case_when(P_value <= 0.001 ~ '***',P_value < 0.01 ~ '**', P_value < 0.05 ~ '*',P_value < 1 ~ 'ns'))

null_full=null_full %>% mutate(color = ifelse(gene %in% c("WT1","WT2","WT3","WT4") , "#89E088", "#88BCE0"))
null_full$gene <- factor(null_full$gene, levels = c( "WT1" ,"WT2" ,"WT3", "WT4","BS1" ,"BS2" ,"BS3", "BS4"))

ggplot(null_full)+geom_violin(aes(gene,enrich,fill=color),lwd=0.9,trim=F) +
	geom_boxplot(width=0.05,aes(gene,enrich),lwd=0.9) +
	#scale_fill_manual(values=c("#53b0db","#53b0db","#53b0db","#6ceb70"))+
	theme_classic(base_size = 15) +
	geom_point(data=true_full,aes(gene,enrich),fill="red",colour="black",pch=21, size=5)+
	theme(legend.position = "none")  +
	geom_text(data=pval,aes(x=File,y=(max(true_full$enrich,null_full$enrich)+0.2),label=sig))  +
	labs(x="Cells",y="Enrichment",title="Genes GRCh37") +
	ggsave("Enrichment/GenomePer/Niek/Niek_Plots/Genes_GRCh37.png")
	

################################################################################
############################## G4s HG38 ######################################
################################################################################
	
null_full=data.frame()
true_full=data.frame()
pval=data.frame()

for (file in list.files("Enrichment/GenomePer/Niek/Niek_G4_hg38//",full.names = T)){
	tmp=read.table(file,header=T)
	true = data.frame(Shift=tmp[1,1],Overlaps=tmp[1,2],gene=strsplit(x=basename(file),split = "[.]")[[1]][1])
	null=tmp[2:nrow(tmp),]
	null$gene=strsplit(x=basename(file),split = "[.]")[[1]][1]
	
	true$enrich = true$Overlaps/mean(null$Overlaps)
	null$enrich = null$Overlaps/mean(null$Overlaps)
	
	null_full=rbind(null_full,null)
	true_full=rbind(true_full,true)
}

pval=data.frame(File=c("BS1","BS2","BS3","BS4","WT1","WT2","WT3","WT4"),P_value=c(0.462,0.374,0.325,0.058,0.289,0.439,0.268,0.352))

pval=pval %>%mutate(sig = case_when(P_value <= 0.001 ~ '***',P_value < 0.01 ~ '**', P_value < 0.05 ~ '*',P_value < 1 ~ 'ns'))

null_full=null_full %>% mutate(color = ifelse(gene %in% c("WT1","WT2","WT3","WT4") , "#89E088", "#88BCE0"))
null_full$gene <- factor(null_full$gene, levels = c( "WT1" ,"WT2" ,"WT3", "WT4","BS1" ,"BS2" ,"BS3", "BS4"))

ggplot(null_full)+geom_violin(aes(gene,enrich,fill=color),lwd=0.9,trim=F) +
	geom_boxplot(width=0.05,aes(gene,enrich),lwd=0.9) +
	#scale_fill_manual(values=c("#53b0db","#53b0db","#53b0db","#6ceb70"))+
	theme_classic(base_size = 15) +
	geom_point(data=true_full,aes(gene,enrich),fill="red",colour="black",pch=21, size=5)+
	theme(legend.position = "none")  +
	geom_text(data=pval,aes(x=File,y=(max(true_full$enrich,null_full$enrich)+0.2),label=sig))  +
	labs(x="Cells",y="Enrichment",title="G4s GRCh38") +
	ggsave("Enrichment/GenomePer/Niek/Niek_Plots/G4s_GRCh38.png")


################################################################################
############################## Genes HG38 ######################################
################################################################################

null_full=data.frame()
true_full=data.frame()
pval=data.frame()

for (file in list.files("Enrichment/GenomePer/Niek/Niek_Genes_hg38///",full.names = T)){
	tmp=read.table(file,header=T)
	true = data.frame(Shift=tmp[1,1],Overlaps=tmp[1,2],gene=strsplit(x=basename(file),split = "[.]")[[1]][1])
	null=tmp[2:nrow(tmp),]
	null$gene=strsplit(x=basename(file),split = "[.]")[[1]][1]
	
	true$enrich = true$Overlaps/mean(null$Overlaps)
	null$enrich = null$Overlaps/mean(null$Overlaps)
	
	null_full=rbind(null_full,null)
	true_full=rbind(true_full,true)
}

pval=data.frame(File=c("BS1","BS2","BS3","BS4","WT1","WT2","WT3","WT4"),
				P_value=c(0.087,0.001,0.038,0.088,0.163,0.425,0.371,0.395))

pval=pval %>%
	mutate(sig = case_when(P_value <= 0.001 ~ '***',
						   P_value < 0.01 ~ '**',
						   P_value < 0.05 ~ '*',
						   P_value < 1 ~ 'ns'))

null_full=null_full %>% mutate(color = ifelse(gene %in% c("WT1","WT2","WT3","WT4") , "#89E088", "#88BCE0"))

null_full$gene <- factor(null_full$gene, levels = c( "WT1" ,"WT2" ,"WT3", "WT4","BS1" ,"BS2" ,"BS3", "BS4"))

ggplot(null_full)+geom_violin(aes(gene,enrich,fill=color),lwd=0.9,trim=F) +
	geom_boxplot(width=0.05,aes(gene,enrich),lwd=0.9) +
	theme_classic(base_size = 15) +
	geom_point(data=true_full,aes(gene,enrich),fill="red",colour="black",pch=21, size=5)+
	theme(legend.position = "none")  +
	geom_text(data=pval,aes(x=File,y=(max(true_full$enrich,null_full$enrich)+0.2),label=sig))  +
	labs(x="Cells",y="Enrichment",title="Genes GRCh38") +
	ggsave("Enrichment/GenomePer/Niek/Niek_Plots/Genes_GRCh38.png")
	
	
