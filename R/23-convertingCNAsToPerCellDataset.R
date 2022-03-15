###################################################
# Makinig CNA per  cell dataset #
################################################### 
libraries=as.character(read.table("Metrics/05-pe.se.metrics.good.quality.ploidy.txt",header=T)$file)
dups=read.table("CNAs/Data/04-somatic_duplications.txt",header=T) %>% select(c(seqnames,start,end,width,gene,file))
dels=read.table("CNAs/03-somatic_deletions.txt",header = T) %>% select(c(seqnames,start,end,width,gene,file))
dups_w=read.table("CNAs/Data/04-somatic_duplications_wt.txt",header=T) %>% select(c(seqnames,start,end,width,gene,file))
dels_w=read.table("CNAs/Data/03-somatic_deletions_wt.txt",header = T) %>% select(c(seqnames,start,end,width,gene,file))
dups=rbind(dups,dups_w)

################################################### 
dups = as.data.frame(dups %>% group_by(file) %>% dplyr::summarize(n()))
colnames(dups)=c("library","cnas")
no_dups= setdiff(libraries,as.character(dups$library))
no_dups_df = data.frame(library=no_dups,cnas=0)
dups = rbind(dups,no_dups_df)
dups$type = "Duplications"

################################################### 
dels=rbind(dels,dels_w)
dels = as.data.frame(dels %>% group_by(file) %>% dplyr::summarize(n()))
colnames(dels)=c("library","cnas")
no_dels= setdiff(libraries,as.character(dels$library))
no_dels_df = data.frame(library=no_dels,cnas=0)
dels = rbind(dels,no_dels_df)
dels$type = "Deletions"

bind = rbind(dels,dups)

################################################### 

bind$library <- bind$library
suppressWarnings(bind <- bind %>% separate(library, c("a","b","c","d","e","f"), "[_-]+"))
bind$gene <- "gene"
for (row in 1:nrow(bind)){
  
  for (letter in c("a","b","c","d","e","f")){
    if (is.na(bind[row,letter])!=T){
      if (bind[row,letter]=="WT" | bind[row,letter]=="wt"){
        bind[row,"gene"]="WT"
      }
      else if (bind[row,letter]=="blm" | bind[row,letter]=="BLM"  | bind[row,letter]=="blm1" ) {
        bind[row,"gene"]="BLM"
      }
      else if (bind[row,letter]=="wrn"  ) {
        if (bind[row,"a"]=="recql5"){
          bind[row,"gene"]="WRN/RECQL5"
        }
        else{
          bind[row,"gene"]="WRN"
        }
      }
      
      else if (bind[row,letter]=="RECQL5" | bind[row,letter]=="recql5" | bind[row,letter]=="RECQ5" | bind[row,letter]=="recq5" ) {
        if (bind[row,"gene"]=="BLM"){
          bind[row,"gene"]="BLM/RECQL5"
        }
        else{
          bind[row,"gene"]="RECQL5"
        }
      }
      else if (bind[row,letter]=="blmrecq5"){
        bind[row,"gene"]="BLM/RECQL5"
      }				
      else if (bind[row,letter]=="RECQL1"| bind[row,letter]=="recql1"){
        bind[row,"gene"]="RECQL1"
      }
      else if ( bind[row,letter]=="fucciwt"| bind[row,letter]=="kbm7" | bind[row,letter]=="wtfucci" ){
        bind[row,"gene"]="WT"
      }
      else if (bind[row,letter]=="RTEL"| bind[row,letter]=="rtel"){
        bind[row,"gene"]="RTEL1"
      }
    }
  }
  
}
bind <- select(bind,-c(a,b,c,d,e,f))
bind$gene=as.factor(bind$gene)

################################################### 

write.table(bind,"CNAs/Data/05-cnas_per_cell_datset_somatic.txt",quote = F,row.names = F,col.names = T,sep="\t")




