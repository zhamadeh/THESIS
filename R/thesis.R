metrics = read.table("Desktop/Metrics.txt",header=T)

metrics$Date=as.Date(metrics$Date,"%d-%m-%Y")
levels(as.factor(metrics$Date))

length(levels(as.factor(metrics$Library)))
nrow(metrics)

n_occur <- data.frame(table(metrics$Library))

n_occur[n_occur$Freq > 1,]

rep=metrics[metrics$Library %in% n_occur$Var1[n_occur$Freq > 1],]

repeated=metrics[metrics$Date=="2020-08-28",]
nrow(repeated)
unique=unique(repeated)

libraries=read.table("Desktop/list.of.thesiis.bam.txt")

setdiff(libraries$V1,paste0(metrics$Library,".trimmed.mdup.bam"))

		