pair = read.table("/Users/zeidh/Desktop/THESIS/Metrics/qualitiyMetricsAllLibraries_pe_se.txt",header = T)

#######################################################
####     Workflow Part 5: Plot Pairs      ###
#######################################################

library(stringr)
library(dplyr)
library("ggplot2")
# Load ggplot2 package

merge =  pair[sample(nrow(pair), 200), ]
	
	#merge=filter(pair, quality_manual != "E")
	merge$quality<-as.factor(droplevels(merge$quality_manual))
	colnames(merge)=c("file" ,"quality_manual","coverage","spikiness", "evenness" ,"evenness.med",
					  "end" ,"background" ,"reads.per.mb.bpr" ,"coverage.bpr" ,"gene" ,"quality"   )
	merge$quality
	my_cols=c("#c92626", "#32a852","#c98d26")
	
	panel.cor <- function(x, y){
		usr <- par("usr"); on.exit(par(usr))
		par(usr = c(0, 1, 0, 1))
		r <- round(cor(x, y), digits=2)
		txt <- paste0("R = ", r)
		cex.cor <- 0.8/strwidth(txt)
		text(0.5, 0.5, txt, cex = cex.cor * r)
	}
	# Customize upper panel
	upper.panel<-function(x, y){
		points(x,y, pch = 21, bg = my_cols[merge$quality],cex=1.4)
	}
	
	panel.hist <- function(x, ...)
	{
		usr <- par("usr"); on.exit(par(usr))
		par(usr = c(usr[1:2], 0, 1.5) )
		h <- hist(x, plot = FALSE)
		breaks <- h$breaks; nB <- length(breaks)
		y <- h$counts; y <- y/max(y)
		rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
	}
	
	# Create the plots
	pairs(merge[, c("coverage","spikiness" , "evenness" , "background" )],
		  lower.panel = panel.cor,
		  upper.panel = upper.panel,
		  diag.panel = panel.hist)
	jpeg(filename = "Output/pairs.png")
	
