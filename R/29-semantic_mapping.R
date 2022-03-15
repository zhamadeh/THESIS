library(GOSemSim)
library(pheatmap) #optional
library(org.Hs.eg.db)

Terms=GOres$description[1:50]
#Make a self-by-self matrix of your GO terms
semsim <- mgoSim(Terms, #character vector of your GO terms
                 Terms, #same vector again. You are comparing the terms to themselves
                 semData=godata('org.Hs.eg.db', ont="BP"), #Human
                 measure="Wang", #distance metric
                 combine=NULL)

#Visualize it with a heatmap (any heatmap function will work, I use pheatmap):
clusters = 5 #The number of semantic clusters to draw. You will want to play around with this number and see what number looks the most correct
pheatmap(semsim, 
         clustering_distance_rows=dist(semsim),
         clustering_distance_cols=dist(semsim),
         border_color = "black",
         treeheight_row = 0,
         cutree_cols = clusters,
         cutree_rows = clusters,
         show_rownames=T,
         color = brewer.pal(9, "BuGn"), #whatever color you want the heatmap to be
         fontsize = 8,
         main="Semantic Clusters")#or whatever

#You can visualize the semsim data however you want really. Dendrogram?? Once you have your clusters, you can use the GOres dataframe to figure out how many of your important genes are annotated in each semantic cluster at least once by:
# for each cluster, extract the associated rows of GOres
# use toString() and unlist() and unique() on the column of GOres containing annotated genes to make a character vector of all the genes in that cluster
# length(important_genes(subsample%in%important_genes))/length(important_genes) to get a percentage, length(important_genes(subsample%in%important_genes)) to get a total count.