ash = read.table("all_bam_prediction.tsv",header = T)
met = read.table("Metrics/qualitiyMetricsAllLibraries_pe_se.txt",header = T) %>% dplyr::select(c(file,quality_manual))
merge = (merge(ash,met,by.x="cell",by.y="file"))
good=filter(merge,merge$quality_manual=="G")
good=levels(droplevels(good$cell))
bad = filter(merge,merge$quality_manual!="G")
bad=levels(droplevels(bad$cell))

df = data.frame(threshold=c(1:100),TP=rep(NA,100),FP=rep(NA,100),FN=rep(NA,100),TN=rep(NA,100))
for (i in 1:100){
	threshold = i/100
	
	tmp = filter(merge,probability > threshold)
	tmp=levels(droplevels(tmp$cell))
	tmp2 = filter(merge,probability <= threshold)
	tmp2=levels(droplevels(tmp2$cell))

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





library(ROCR)

df=dplyr::select(merge,c(probability,quality_manual))
df$quality_manual=gsub(x=df$quality_manual,pattern = "G",replacement = 1)
df$quality_manual=gsub(x=df$quality_manual,pattern = "E",replacement = 0)
df$quality_manual=gsub(x=df$quality_manual,pattern = "P",replacement = 0)
pred <- prediction(df$probability, df$quality_manual)
df_pred = data.frame(TP=pred@tp[[1]],FP=pred@fp[[1]],TN=pred@tn[[1]],FN=pred@fn[[1]])
df_pred$sensitivity=(df_pred$TP/(df_pred$TP+df_pred$FN))
df_pred$specificity=1-(df_pred$TN/(df_pred$TN+df_pred$FP))
df_pred$accuracy = (df_pred$TP+df_pred$TN)/(df_pred$TP+df_pred$FP+df_pred$TN+df_pred$FN)
max(df_pred$accuracy)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE)
results(perf)

df = data.frame(TPR=perf@y.values[[1]], FPR=perf@x.values[[1]])
ggplot(df)+geom_line(aes(FPR,TPR))








metrics=read.table("/Users/zeidh//Desktop/THESIS/Metrics/qualitiyMetricsAllLibraries_pe_se.txt",header = T)


metrics$quality_manual <- droplevels.factor(metrics$quality_manual)
metrics %>% group_by(quality_manual)%>% summarize((n()))
metrics$quality_manual = gsub(x = metrics$quality_manual,pattern = "E",replacement = "P")

metrics=dplyr::select(metrics,-c(file,gene,end))
################## Data partitioning ##################
# create a list of 80% of the rows in the original dataset we can use for training


validation_index <- createDataPartition(metrics$quality_manual, p=0.80, list=FALSE)
testData <- metrics[-validation_index,] # select 20% of the data for validation
trainData <- metrics[validation_index,] # use the remaining 80% of data to training and testing the models

fit.rf <- train(quality_manual~., data=trainData, method=, metric=metric, trControl=control)

testData1=testData
testData2=dplyr::select(testData,-quality_manual)
predictions <- predict(fit.rf, na.omit(testData2),type = "prob")
head(predictions)
rocCurve.bagg <- roc(testData1$quality_manual,predictions[,"G"])
#plot the ROC curve
rocCurve.bagg$sensitivities
rocCurve.bagg$specificities
plot(rocCurve.bagg,col=c(6))


df1 = data.frame(TPR=perf@y.values[[1]], FPR=perf@x.values[[1]])
df1$type="ASHLEYs (Gros et al., 2021)"
df2 = data.frame(TPR=rocCurve.bagg$sensitivities, FPR=1-rocCurve.bagg$specificities)
df2$type="New Method"
df=rbind(df1,df2)
ggplot(df)+geom_line(aes(FPR,TPR,color=type),size=1.3)+
	xlim(0,1)+ylim(0,1)+
	labs(x="False Positive Rate (1 - Specificity)",y="True Positive Rate (Sensitivity)")+
	theme_classic()+
	theme(text = element_text(size = 16),aspect.ratio = 1,legend.position = c(0.8, 0.3))+
	guides(color=guide_legend(title="Classifier"))+
	annotate(geom="text", x=0.3, y=0.96, label="Accuracy = 96.6%",
			   color="#00BFC4")+
	annotate(geom="text", x=0.74, y=0.84, label="Accuracy = 83.4%",
			 color="#F8766D")+
	ggsave("Output/ROC_Ashleys_vs_newMethod.png")

df$accuracy = (df$TPR+df$FPR)/2

df %>% group_by(type)  %>% summarize(max())
