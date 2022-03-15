######################################################################
####   Plotting Karyogram   ####
######################################################################

################ Step 1: load genomic range data ####################

# Needs to be in dataframe format, can be points or ranges
SNP=read.table('CNAs/somatic_cnas.txt',header=T) %>% select(c(seqnames,start,end,type,gene))
#cnvPerCell <- read.table(paste0("CNVs/FUCCI/30Mb_cnvPerCell.txt"),header=T)
#SNP=filter(cnvPerCell,type=='gain') %>% select(c(seqnames,start,end))

# Give bed file column names for first three rows
colnames(SNP)=c("chr","start","end", "type"   ,  "gene"  )

# Convert chrom to factor
SNP$chr<-as.factor(SNP$chr)

# Reorder levels for printing 
SNP$chr = factor(SNP$chr,levels=c("chr1" ,"chr2" ,"chr3" , "chr4"  ,"chr5" , "chr6"  ,"chr7"  ,"chr8",  "chr9" ,
								  "chr10", "chr11" ,"chr12", "chr13" ,"chr14" ,"chr15" ,"chr16" ,"chr17", "chr18" ,"chr19" , 
								  "chr20", "chr21","chr22" , "chrX" , "chrY" ))

# Divide by 1Mb to simplify axis in plotting
SNP$start=SNP$start/1000000
SNP$end = SNP$end/1000000


################ Step 2: load Genome data (lengths) ####################

# Load chromosome lengths
data = read.table("Anxilliary/hg38_chrLengths.txt") 

# Make sure they have "chr" nomenclature
data$V1<-paste0("chr",data$V1)

# Get rid of 'start' column you only need the lengths
data=select(data,-c(V2))

# Rename columns to chromosome and size
colnames(data)=c("chr","size")

# Change to factor
data$chr<-as.factor(data$chr)

# Reorder levels for plotting
data$chr = factor(data$chr,levels=c("chr1" ,"chr2" ,"chr3" , "chr4"  ,"chr5" , "chr6"  ,"chr7"  ,"chr8",  "chr9" ,
									"chr10", "chr11" ,"chr12", "chr13" ,"chr14" ,"chr15" ,"chr16" ,"chr17", "chr18" ,"chr19" , 
									"chr20", "chr21","chr22" , "chrX" , "chrY" ))

# Divide by 1Mb to clean up axis
data$size=data$size/1000000
data[data$chr=="chr9",]$size=data[data$chr=="chr9",]$size/2


################ Step 3: load secondary dataset (eg. CFSs) ####################

# Load data
cfs = read.table("Anxilliary/Fragile_sites/hg38_cfs.bed")

# Change lowercase 'x'
cfs$V1=gsub("chrx","chrX",cfs$V1)

# Divide by 1Mb for axis
cfs$V3=cfs$V3/1000000
cfs$V2=cfs$V2/1000000

# Change to factor and reorder levels
cfs$V1<-as.factor(cfs$V1)
cfs$V1 = factor(cfs$V1,levels=c("chr1" ,"chr2" ,"chr3" , "chr4"  ,"chr5" , "chr6"  ,"chr7"  ,"chr8",  "chr9" ,
								"chr10", "chr11" ,"chr12", "chr13" ,"chr14" ,"chr15" ,"chr16" ,"chr17", "chr18" ,"chr19" , 
								"chr20", "chr21","chr22" , "chrX" , "chrY" ))



################ Step 4: Plotting  ####################

ggplot() +
	geom_segment(data = data,
				 aes(x = chr, xend = chr, y = 0, yend = size),
				 lineend = "round", color = "lightgrey", size = 5) +
	geom_rect(data = SNP,
			  aes(xmin = as.integer(chr) - 0.25, xmax = as.integer(chr) + 0.25,
			  	ymin = start, ymax = end)
			  ,fill="#1d705d",size = 0.25,alpha=0.2) +
	#geom_rect(data = cfs, aes(xmin = as.integer(V1) - 0.4, xmax = as.integer(V1) - 0.3,
	#						  ymin = V2, ymax = V3) ,fill="black",
	#		  size = 0.25) +
	theme_classic() +
	theme(text = element_text(size=18),axis.line=element_blank())+
	labs(x="",y="CHROMOSOME POSITION (Mb)")#+
	#ggsave("OUTPUT/Karyograms/karyogram_curatedDeletions.png")




######################################################################
####   Plotting Karyogram   ####
######################################################################

################ Step 1: load genomic range data ####################

# Needs to be in dataframe format, can be points or ranges
SNP1=read.table("CNAs/somatic_duplications.txt",header=T) %>% select(c(seqnames,start,end,gene))
SNP2=read.table("CNAs/somatic_deletions.txt",header = T) %>% select(c(seqnames,start,end,gene))
#cnvPerCell <- read.table(paste0("CNVs/FUCCI/30Mb_cnvPerCell.txt"),header=T)
#SNP=filter(cnvPerCell,type=='gain') %>% select(c(seqnames,start,end))

# Give bed file column names for first three rows
colnames(SNP1)=c("chr","start","end","gene")
colnames(SNP2)=c("chr","start","end","gene")

# Convert chrom to factor
SNP1$chr<-as.factor(SNP1$chr)
SNP2$chr<-as.factor(SNP2$chr)

# Reorder levels for printing 
SNP1$chr = factor(SNP1$chr,levels=c("chr1" ,"chr2" ,"chr3" , "chr4"  ,"chr5" , "chr6"  ,"chr7"  ,"chr8",  "chr9" ,
									"chr10", "chr11" ,"chr12", "chr13" ,"chr14" ,"chr15" ,"chr16" ,"chr17", "chr18" ,"chr19" , 
									"chr20", "chr21","chr22" , "chrX" , "chrY" ))
# Reorder levels for printing 
SNP2$chr = factor(SNP2$chr,levels=c("chr1" ,"chr2" ,"chr3" , "chr4"  ,"chr5" , "chr6"  ,"chr7"  ,"chr8",  "chr9" ,
									"chr10", "chr11" ,"chr12", "chr13" ,"chr14" ,"chr15" ,"chr16" ,"chr17", "chr18" ,"chr19" , 
									"chr20", "chr21","chr22" , "chrX" , "chrY" ))

# Divide by 1Mb to simplify axis in plotting
SNP1$start=SNP1$start/1000000
SNP1$end = SNP1$end/1000000
SNP2$start=SNP2$start/1000000
SNP2$end = SNP2$end/1000000


################ Step 2: load Genome data (lengths) ####################

# Load chromosome lengths
data = read.table("Anxilliary/hg38_chrLengths.txt") 

# Make sure they have "chr" nomenclature
data$V1<-paste0("chr",data$V1)

# Get rid of 'start' column you only need the lengths
data=select(data,-c(V2))

# Rename columns to chromosome and size
colnames(data)=c("chr","size")

# Change to factor
data$chr<-as.factor(data$chr)

# Reorder levels for plotting
data$chr = factor(data$chr,levels=c("chr1" ,"chr2" ,"chr3" , "chr4"  ,"chr5" , "chr6"  ,"chr7"  ,"chr8",  "chr9" ,
									"chr10", "chr11" ,"chr12", "chr13" ,"chr14" ,"chr15" ,"chr16" ,"chr17", "chr18" ,"chr19" , 
									"chr20", "chr21","chr22" , "chrX" , "chrY" ))

# Divide by 1Mb to clean up axis
data$size=data$size/1000000
data[data$chr=="chr9",]$size=data[data$chr=="chr9",]$size/2



################ Step 3: load secondary dataset (eg. CFSs) ####################

# Load data
cfs = read.table("Enrichment/hg38_cfs.bed")
cfs$V1=paste0("chr",cfs$V1)
# Change lowercase 'x'
cfs$V1=gsub("chrx","chrX",cfs$V1)

# Divide by 1Mb for axis
cfs$V3=cfs$V3/1000000
cfs$V2=cfs$V2/1000000

# Change to factor and reorder levels
cfs$V1<-as.factor(cfs$V1)
cfs$V1 = factor(cfs$V1,levels=c("chr1" ,"chr2" ,"chr3" , "chr4"  ,"chr5" , "chr6"  ,"chr7"  ,"chr8",  "chr9" ,
								"chr10", "chr11" ,"chr12", "chr13" ,"chr14" ,"chr15" ,"chr16" ,"chr17", "chr18" ,"chr19" , 
								"chr20", "chr21","chr22" , "chrX" , "chrY" ))



################ Step 4: Plotting  ####################
#wt=data.frame(chr=c("chr8","chr15"),start=c(0.000001,60.8),end=c(139.88456,89.4),gene=c("WT","WT"))
#SNP1=rbind(SNP1,wt)


ggplot() +
	geom_segment(data = data,
				 aes(x = chr, xend = chr, y = 0, yend = size),
				 lineend = "round", color = "lightgrey", size = 5) +
	geom_rect(data = SNP1,
			  aes(xmin = as.integer(chr) - 0.25, xmax = as.integer(chr),
			  	ymin = start, ymax = end)
			  ,fill="#036bfc",size = 0.25,alpha=0.6) +
	geom_rect(data = SNP2,
			  aes(xmin = as.integer(chr) , xmax = as.integer(chr) + 0.25,
			  	ymin = start, ymax = end)
			  ,fill="#fc1c03",size = 0.25,alpha=0.6) +
	#geom_rect(data = cfs, aes(xmin = as.integer(V1) - 0.34, xmax = as.integer(V1) - 0.27,
	#                         ymin = V2, ymax = V3) ,fill="black",
	#          size = 0.25) +
	theme_classic() +
	theme(text = element_text(size=21),axis.line=element_blank(),
		  axis.text.x=element_blank(),
		  axis.ticks.x=element_blank())+
	labs(x="Chromosomes",y="Chromosome position (Mb)")+
	facet_wrap(~gene,ncol = 1)+
	ggsave("Output/Del_Dup_overlay_ByGene_karyogram.png")
