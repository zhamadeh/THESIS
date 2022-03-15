#################################################################################
#################################################################################

require("caret")
require("Boruta")
require("RANN")
library(tidyverse)
library(caret)
library(Boruta)
library(tidyr)
library(RANN)
library(ROCR)
library(pROC)

#################################################################################
############################## Training model ##############################
#################################################################################
metrics=read.table("/Users/zeidh//Desktop/THESIS/Metrics/qualitiyMetricsAllLibraries_pe_se.txt",header = T)
metrics=select(metrics,-c(evenness.med,coverage.bpr,ploidy))
metrics$quality_manual <- droplevels.factor(metrics$quality_manual)
metrics %>% group_by(quality_manual)%>% summarize((n()))
metrics$quality_manual = gsub(x = metrics$quality_manual,pattern = "E",replacement = "P")
metrics$quality_manual=as.factor(metrics$quality_manual)
metrics=dplyr::select(metrics,-c(file,gene,end))
# Data partitioning
# create a list of 80% of the rows in the original dataset we can use for training
control <- trainControl(method="cv", number=10,classProbs = T)
metric <- "Accuracy"

validation_index <- createDataPartition(metrics$quality_manual, p=0.80, list=FALSE)
testData <- metrics[-validation_index,] # select 20% of the data for validation
trainData <- metrics[validation_index,] # use the remaining 80% of data to training and testing the models

fit.rf <- train(quality_manual~., data=trainData, method="rf", metric=metric, trControl=control)
predictions <- predict(fit.rf, na.omit(testData),type = "prob")
#confusionMatrix(predictions, testData$quality_manual)
rocCurve.bagg_zh_lansdorp <- roc(testData$quality_manual,predictions[,"G"])

#################################################################################
############################## Testing ploidy ##############################
#################################################################################
metrics=read.table("/Users/zeidh//Desktop/THESIS/Metrics/qualitiyMetricsAllLibraries_pe_se.txt",header = T)
metrics=filter(metrics,ploidy==2)
metrics=select(metrics,-c(evenness.med,coverage.bpr,ploidy,file,gene,end))
metrics$quality_manual <- droplevels.factor(metrics$quality_manual)
metrics$quality_manual = gsub(x = metrics$quality_manual,pattern = "E",replacement = "P")
metrics$quality_manual=as.factor(metrics$quality_manual)
control <- trainControl(method="cv", number=10,classProbs = T)
metric <- "Accuracy"

validation_index <- createDataPartition(metrics$quality_manual, p=0.80, list=FALSE)
testData <- metrics[-validation_index,] # select 20% of the data for validation
trainData <- metrics[validation_index,] # use the remaining 80% of data to training and testing the models

fit.rf <- train(quality_manual~., data=trainData, method="rf", metric=metric, trControl=control)
predictions <- predict(fit.rf, na.omit(testData),type = "prob")
rocCurve.bagg_zh_lansdorp_ploidy2 <- roc(testData$quality_manual,predictions[,"G"])

################################# n = 1################################################

metrics=read.table("/Users/zeidh//Desktop/THESIS/Metrics/qualitiyMetricsAllLibraries_pe_se.txt",header = T)
metrics=filter(metrics,ploidy==1)
metrics=select(metrics,-c(evenness.med,coverage.bpr,ploidy,file,gene,end))
metrics$quality_manual <- droplevels.factor(metrics$quality_manual)
metrics$quality_manual = gsub(x = metrics$quality_manual,pattern = "E",replacement = "P")
metrics$quality_manual=as.factor(metrics$quality_manual)
control <- trainControl(method="cv", number=10,classProbs = T)
metric <- "Accuracy"

validation_index <- createDataPartition(metrics$quality_manual, p=0.80, list=FALSE)
testData <- metrics[-validation_index,] # select 20% of the data for validation
trainData <- metrics[validation_index,] # use the remaining 80% of data to training and testing the models

fit.rf <- train(quality_manual~., data=trainData, method="rf", metric=metric, trControl=control)
predictions <- predict(fit.rf, na.omit(testData),type = "prob")
rocCurve.bagg_zh_lansdorp_ploidy1 <- roc(testData$quality_manual,predictions[,"G"])

#################################################################################
############################## Making predictioins ##############################
#################################################################################

ash = read.table("Classifier/NA878 test sets/na878-porubsky/na878-porubsky-classifier-test-metrics.txt",header = T) %>% dplyr::select(-c(evenness.med,coverage.bpr))
porubsky=dplyr::select(ash,-c(file,ASHLEYS))
predictions <- predict(fit.rf, porubsky,type = "prob")
rocCurve.bagg_zh_porubsky <- roc(ash$quality_manual,predictions[,"G"])

#################################################################################
############################## Making predictioins ##############################
#################################################################################

df=dplyr::select(ash,c(ASHLEYS,quality_manual))
df$quality_manual=gsub(x=df$quality_manual,pattern = "G",replacement = 1)
df$quality_manual=gsub(x=df$quality_manual,pattern = "E",replacement = 0)
df$quality_manual=gsub(x=df$quality_manual,pattern = "P",replacement = 0)
pred <- prediction(df$ASHLEYS, df$quality_manual)
rocCurve.bagg_ash_porubsky <- roc(ash$quality_manual,pred@predictions[[1]])

#run on lansdorp data
ash = read.table("Classifier/FUCCI test sets/all_bam_prediction.tsv",header = T)
met = read.table("Metrics/qualitiyMetricsAllLibraries_pe_se.txt",header = T) %>% dplyr::select(c(file,quality_manual))
merge = (merge(ash,met,by.x="cell",by.y="file"))
df=dplyr::select(merge,c(probability,quality_manual))
df$quality_manual=gsub(x=df$quality_manual,pattern = "G",replacement = 1)
df$quality_manual=gsub(x=df$quality_manual,pattern = "E",replacement = 0)
df$quality_manual=gsub(x=df$quality_manual,pattern = "P",replacement = 0)
pred <- prediction(df$probability, df$quality_manual)
perf <- performance(pred,"tpr","fpr")
#rocCurve.bagg_ash_lansdorp = data.frame(TPR=perf@y.values[[1]], FPR=perf@x.values[[1]])
rocCurve.bagg_ash_lansdorp <- roc(df$quality_manual,pred@predictions[[1]])

pred@predictions

#################################################################################
############################## Making predictioins ##############################
#################################################################################

df1 = data.frame(TPR=rocCurve.bagg_ash_porubsky$sensitivities, FPR=1-rocCurve.bagg_ash_porubsky$specificities)
df1$type="ASHLEYs (Gros et al., 2021)"
df1$data="Porubsky (2016)"
df2 = data.frame(TPR=perf@y.values[[1]], FPR=perf@x.values[[1]])
df2$type="ASHLEYs (Gros et al., 2021)"
df2$data="Lansdorp Lab (2020)"

df3 = data.frame(TPR=rocCurve.bagg_zh_porubsky$sensitivities, FPR=1-rocCurve.bagg_zh_porubsky$specificities)
df3$type="New Method"
df3$data="Porubsky (2016)"
df4 = data.frame(TPR=rocCurve.bagg_zh_lansdorp$sensitivities, FPR=1-rocCurve.bagg_zh_lansdorp$specificities)
df4$type="New Method"
df4$data="Lansdorp Lab (2020)"


df=rbind(df1,df2,df3,df4)
ggplot(df)+geom_line(aes(FPR,TPR,color=type,linetype=data),size=1.3)+
	xlim(0,1)+ylim(0,1.05)+
	labs(x="False Positive Rate (1 - Specificity)",y="True Positive Rate (Sensitivity)")+
	theme_classic()+
	theme(text = element_text(size = 16),aspect.ratio = 1,legend.position = c(0.6, 0.3))+
	guides(color=guide_legend(title="Classifier"),linetype=guide_legend(title="Test Data"))+
	annotate(geom="text", x=0.35, y=0.9, label=paste0("AUC = ",round(rocCurve.bagg_ash_porubsky$auc[[1]],digits = 2)),
			 color="#F8766D")+
	annotate(geom="text", x=0.78, y=0.84, label=paste0("AUC = ",round(rocCurve.bagg_ash_lansdorp$auc[[1]],digits = 2)),
			 color="#F8766D")+
	annotate(geom="text", x=0.21, y=0.70, label=paste0("AUC = ",round(rocCurve.bagg_zh_porubsky$auc[[1]],digits = 2)),
			 color="#00BFC4")+
	annotate(geom="text", x=0.1, y=1.03, label=paste0("AUC = ",round(rocCurve.bagg_zh_lansdorp$auc[[1]],digits = 2)),
			 color="#00BFC4")+ggsave("Output/ROC_Ashleys_vs_newMethod.png")


#################################################################################
#################################################################################

df1 = data.frame(TPR=rocCurve.bagg_zh_lansdorp_ploidy1$sensitivities, FPR=1-rocCurve.bagg_zh_lansdorp_ploidy1$specificities)
df1$type=paste0("Haploid (AUC = ",round(rocCurve.bagg_zh_lansdorp_ploidy1$auc[[1]],digits = 2),")")
df2 = data.frame(TPR=rocCurve.bagg_zh_lansdorp_ploidy2$sensitivities, FPR=1-rocCurve.bagg_zh_lansdorp_ploidy2$specificities)
df2$type=paste0("Diploid (AUC = ",round(rocCurve.bagg_zh_lansdorp_ploidy2$auc[[1]],digits = 2),")")
df=rbind(df1,df2)
ggplot(df)+geom_line(aes(FPR,TPR,color=type),size=1.3)+
	xlim(0,1)+ylim(0,1)+
	labs(x="False Positive Rate (1 - Specificity)",y="True Positive Rate (Sensitivity)")+
	theme_classic()+
	theme(text = element_text(size = 16),aspect.ratio = 1,legend.position = c(0.8, 0.3))+
	guides(color=guide_legend(title="Classifier"))
	#annotate(geom="text", x=0.74, y=0.60, label=paste0("AUC = ",round(rocCurve.bagg_zh_lansdorp_ploidy2$auc[[1]],digits = 2)),
			 #color="#F8766D")+ 
	+ggsave("Output/ROC_newMethod_ploidybased.png")


#########################################################################################








merge=ash
good=filter(merge,merge$quality_manual=="G")
good=levels(droplevels(good$file))
bad = filter(merge,merge$quality_manual!="G")
bad=levels(droplevels(bad$file))

df = data.frame(threshold=c(1:100),TP=rep(NA,100),FP=rep(NA,100),FN=rep(NA,100),TN=rep(NA,100))
for (i in 1:100){
	threshold = i/100
	
	tmp = filter(merge,ASHLEYS > threshold)
	tmp=levels(droplevels(tmp$file))
	tmp2 = filter(merge,ASHLEYS <= threshold)
	tmp2=levels(droplevels(tmp2$file))
	
	df[i,2]=length(intersect(tmp,good))
	df[i,3]=length(setdiff(tmp,good))
	df[i,4]=length(setdiff(good,tmp))
	df[i,5]=length(intersect(bad,tmp2))
	
	
} 
df$good=length(good)
df$TPpFN=df$TP+df$FN
df$bad=length(bad)
df$TNpFP=df$TN+df$FP

df$sensitivity=(df$TP/(df$TP+df$FN))
df$specificity=1-(df$TN/(df$TN+df$FP))
df$good=length(good)
df$bad=length(bad)
df


ggplot(df)+geom_point(aes(specificity,sensitivity))+
	geom_smooth(aes(specificity,sensitivity),se=F)+
	xlim(0,1)+ylim(0,1)+
	labs(x="False Positive Rate (1 - Specificity)",y="True Positive Rate (Sensitivity)")+
	theme_classic()+
	theme(text = element_text(size = 15),aspect.ratio = 1)#+
#	ggsave("Precision Recall/prec_recall.png")

