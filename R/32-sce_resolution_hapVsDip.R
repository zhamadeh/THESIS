
########################################################################
################## SCEs from haploid and diploid cells ################
########################################################################

fucci=read.table("SCEs/RESOLUTION/haploidSCEs_Resolution.txt",header=T)
merge=read.table("SCEs/RESOLUTION/diploidSCEs_Resolution.txt",header=T)

fucci$type=paste0("Haploid (n = ",nrow(fucci),")")
merge$type=paste0("Diploid (n = ",nrow(merge),")")

bind=rbind(fucci,merge)


ggplot(bind)+geom_smooth(aes(reads.per.mb.bpr,width,color=type))+
	theme_classic()+
	#geom_text(data=pval,aes(x=gene,y=(max(test$norm)+5),label=sig)) +
	theme(text=element_text(size=15),legend.position = c(0.6,0.7))+
	guides(color=guide_legend(title="Cell Ploidy"))+
	labs(x="Sequencing Effort (Reads/Mb)",y="SCE Breakpoint Resolution")+
	annotate(geom="text", x=200, y=9e+04, label=paste0("p = ",round(p_value,digits = 10)))+
	ggsave("Output/SCEs_hap_dip_resolution.png")


### https://stats.stackexchange.com/questions/435644/is-there-a-method-to-look-for-significant-difference-between-two-linear-regressi

compare.coeff <- function(b1,se1,b2,se2){
	return((b1-b2)/sqrt(se1^2+se2^2))
}
lm1 = lm(width ~ reads.per.mb.bpr,data=subset(bind,bind$type=="haploid"))
lm2 = lm(width ~ reads.per.mb.bpr,data=subset(bind,bind$type=="diploid"))

b1 <- summary(lm1)$coefficients[2,1]
se1 <- summary(lm1)$coefficients[2,2]
b2 <- summary(lm2)$coefficients[2,1]
se2 <- summary(lm2)$coefficients[2,2]

p_value = 2*pnorm(-abs(compare.coeff(b1,se1,b2,se2)))
p_value

summary(lm(width ~ reads.per.mb.bpr + type + reads.per.mb.bpr:type,data=bind))

########################################################################
################## Original datasets of SCEs ################
########################################################################

fucci=read.table("SCEs/RESOLUTION/Initial/fucci_sces_complete.bed",header=T)
fucci=filter(fucci,ploidy==1)
fucci=select(fucci,c(width,reads.per.mb.bpr))
colnames(fucci)
fucci$type="haploid"

write.table(fucci,"SCEs/RESOLUTION/haploidSCEs_Resolution.txt",quote = F,row.names = F,col.names = T,sep = "\t")

na=read.table("SCEs/RESOLUTION/Initial/na878_sces.txt",header=T) %>% select(c(start,end,filenames))
na$filenames=tools::file_path_sans_ext(na$filenames)
na$filenames=tools::file_path_sans_ext(na$filenames)
head(na)
nrow(na)
rpm=read.table("SCEs/RESOLUTION/Initial/na878_reads_perMB.txt",header=T)
nrow(rpm)
merge=merge(na,rpm,by.x="filenames",by.y="file")
merge$width=merge$end-merge$start
merge=select(merge,c(width,RPM))
colnames(merge)=c("width"  ,"reads.per.mb.bpr")
merge$type="diploid"

write.table(merge,"SCEs/RESOLUTION/diploidSCEs_Resolution.txt",quote = F,row.names = F,col.names = T,sep = "\t")
