#Written in R 4.0.3 by Kat Lande

#http://www.webgestalt.org/results/1646855873/#

library(stringr)
library(dplyr)
library(ggplot2)
library(WebGestaltR)
library(RColorBrewer)
#remove.packages("Rcpp")
#install.packages("Rcpp")
#update.packages("Rcpp")
library(Rcpp)
#there might be other packages you need -- I'm not sure? These functions come from a much larger script.

#Make a table of GO terms with ORA, NTA, or GSEA
GOtable <- function(res, library, direction=NULL, organism = "hsapiens", 
                    enrichMethod = "ORA", permutes=1000, neighbors=10, 
                    highlightSeedNum=10, topThr=10, mode="df", genecol=7, lfccol=2, pathway="pathway_KEGG"){
  
  
  
  if(enrichMethod == "GSEA"){
    resSig <- res
    gene_list <- row.names(resSig)
  }else{
    if(mode == "df"){
      if(is.null(direction) == F){
        if(direction == "UP"){
          resSig <- subset(res, log2FoldChange > 0)
        }else if(direction == "DOWN"){
          resSig <- subset(res, log2FoldChange < 0)
        }else{print("ERROR: Direction not recognized")}
      } else { resSig <- res}
      gene_list <- row.names(resSig)
    }else if(mode == "vector"){
      gene_list <- res
    }
  }
  
  
  
  print(paste0(length(gene_list), " Genes to test..."))
  
  if(enrichMethod == "ORA"){
    tab <- WebGestaltR(enrichMethod="ORA",
                       organism=organism,
                       enrichDatabase="geneontology_Biological_Process", #All BP GO terms
                       interestGeneType="genesymbol",
                       referenceGeneType="genesymbol",
                       referenceGene=library,
                       interestGene=gene_list)
    
  }else if(enrichMethod == "GSEA"){
    resSig$gene <- row.names(resSig)
    ranks <- resSig[c(genecol,lfccol)]
    colnames(ranks)[2] <- "rank"
    ranks$rank <- as.numeric(unlist(ranks$rank))
    
    
    tab <- WebGestaltR(enrichMethod="GSEA", organism=organism,
                       enrichDatabase=pathway, interestGene=ranks, #uses KEGG pathways
                       interestGeneType="genesymbol", sigMethod="top", 
                       topThr=topThr, minNum=5, perNum = permutes,
                       saveRawGseaResult=TRUE)
    
  }else if(enrichMethod == "NTA"){
    tab <- WebGestaltR(enrichMethod="NTA", 
                       organism=organism,
                       enrichDatabase="network_PPI_BIOGRID", #uses protein-protein interaction networks
                       interestGene=gene_list,
                       interestGeneType="genesymbol", 
                       sigMethod="top", 
                       topThr=topThr,
                       highlightSeedNum=highlightSeedNum, 
                       neighborNum = neighbors,
                       networkConstructionMethod="Network_Retrieval_Prioritization")}
  
  return(tab)
}
# res: df wherein row.names() == genes, or vector of hit genes
# library: a vector of all genes in the genome (the background library for you GO query)
# direction: "UP" or "DOWN", whether to look at logFC positive or negative genes. You can ignore this if you are not looking at RNAseq data.
# organism: organism input for WebGestaltR, auto="hsapiens"
# enrichMethod: gene ontology method; "ORA"/"GSEA"/"NTA", auto = "ORA" 
# permutes: permutations to run for GSEA, auto=1000.
# neighbors: neighborNum for NTA, auto=10, 
# highlightSeedNum: highlightSeedNum for NTA, auto=10
# topThr: number of terms to pull out for NTA/GSEA, auto=10
# mode: whether res is "df" or "vector"
# genecol,lfccol: columns of the "rank" df corresponding to gene and rank (LFC) variable, auto=7,2, only relevant for GSEA. 


#Make a barplot of top ORA terms
ORA.barplot <- function(ORA, top.terms=30, title="ORA Plot", continuous=T, omitNA = T){
  #add a labels column to our data frame that has the GO term & description
  ORA$labs <- ORA$description#paste(ORA$geneSet,ORA$description, sep=": ")
  #create a new df of only siginificant results
  ORA.Sig <- ORA[ORA$FDR < 0.05,]
  #order this df by FDR
  #print(is.numeric(ORA.Sig$FDR)) #Is numeric
  ORA.Sig <- ORA.Sig[order(ORA.Sig$FDR, decreasing= F),]
  #make the label column into an ordered factor 
  ORA.Sig$labs <- factor(ORA.Sig$labs, levels = ORA.Sig$labs[order(-log(ORA.Sig$FDR))])
  if(omitNA == T){
    dat <- na.omit(ORA[1:top.terms,])
  }else{dat <- ORA[1:top.terms,]}
  #no idea why I have to do this twice ?
  dat$labs <- factor(dat$labs, levels = dat$labs[order(-log(dat$FDR))])
  
  plot <- ggplot(data=dat, aes(x=labs, y=(-log(FDR)),fill=enrichmentRatio)) + 
    geom_bar(stat="identity",col="black") +
    geom_text(aes(label=signif(FDR, digits = 2)), hjust=1.1, color="white", size=3) +
    labs(title=title) +
    theme(axis.text.x = element_text(size=1), plot.title = element_text(hjust = 0.5))+
    coord_flip() +
    theme_minimal() +
    theme(axis.title.y=element_blank())
  if(continuous == T){
    plot <- plot + scale_fill_gradient("Observed/Expected", low="#B0E2AF", high="#2B83BA")
  }
  
  
  return(plot)
}
# ORA: ORA results from GOtable()
# top.terms: terms to include in plot, auto=30
# title: plot title, auto="ORA Plot"
# continuous: whether the color variable is continuous, auto=T
# omitNA: T/F, whether to exclude rows containing NA values, auto=T


# Use example for Zeid:
Important_Genes=as.character(read.table("Enrichment/Genes/Initial/essentialKBM7genes.txt",header=T)$Gene)
Background_Genes=as.character(levels(droplevels(as.factor(genes$geneName))))
#1 - find GO terms
GOres <- GOtable(Important_Genes, #All your "important genes," as a character vector.
        Background_Genes, #All genes in the genome. This acts as the "background library" to query against for enrichment. Input is a character vector.
        organism = "hsapiens", #I assume this is human data? change if not
        enrichMethod = "ORA", #ORA is the simplest type of GO. I don't think you can do the other types based on the way to collected your data.
        pathway="geneontology_Biological_Process_noRedundant",#BP terms are the most informative IMO. WebGestaltR has other GO libraries you can query if the results aren't what you expected. Check their documentation for lists.
        mode="vector") 

#2 - plot terms
ORA.barplot(GOres, #output from GOtable()
            top.terms=50, #the number of terms you want to show up on the plot
            title="ORA Plot") #Whatever you want the title to be

write.table(GOres,"GO_analysis.txt",quote = F,row.names = F,col.names = T,sep = "\t")



