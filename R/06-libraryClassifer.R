################### Required packages ##################
#install.packages("caret", dependencies=c("Depends", "Suggests"))
#install.packages('Boruta')
#install.packages('RANN')
require("caret")
require("Boruta")
require("RANN")
library(tidyverse)
library(caret)
library(Boruta)
library(tidyr)
library(RANN)

################### Assembly of training set ##################

metrics=read.table("/Users/zeidh//Desktop/THESIS/Metrics/qualitiyMetricsAllLibraries_pe_se.txt",header = T)


metrics$quality_manual <- droplevels.factor(metrics$quality_manual)
metrics %>% group_by(quality_manual)%>% summarize((n()))
metrics$quality_manual = gsub(x = metrics$quality_manual,pattern = "E",replacement = "P")

metrics=select(metrics,-c(file,gene,end))
################## Data partitioning ##################
# create a list of 80% of the rows in the original dataset we can use for training


validation_index <- createDataPartition(metrics$quality_manual, p=0.80, list=FALSE)
testData <- metrics[-validation_index,] # select 20% of the data for validation
trainData <- metrics[validation_index,] # use the remaining 80% of data to training and testing the models

#dataset=na.exclude(dataset)


################## Pre processing ###################

#Note imputation made performance on validation set worse even though there are no NAs
#Meaning norrmalizatiton was having a negative effect

#preProcess_missingdata_model <- preProcess(trainData, method='knnImpute')
#trainData <- predict(preProcess_missingdata_model, newdata = trainData)  #  !anyNA(trainData)
#preProcess_missingdata_model <- preProcess(testData, method='knnImpute')
#testData <- predict(preProcess_missingdata_model, newdata = testData)

#Note if we had any categorical variables we would need to one hot encode them
#But in this case we don't so our dataset is complete

################### Preliminary plotting ##################
#just  to see if any variables are obviously important

#x <- trainData[,c("Reads","complexity","Background","spikiness","entropy","coverage")] #metrics
#y <- trainData[,5] #classification
#scales <- list(x=list(relation="free"), y=list(relation="free"))
#featurePlot(x=x, y=y, plot="box",scales =scales)
#featurePlot(x=x, y=y, plot="density",scales =scales,strip=strip.custom(par.strip.text=list(cex=.7)))
#featurePlot(x = metrics[, c("Background","Reads","complexity","spikiness","entropy")],y = metrics$quality_manual,plot = "pairs",auto.key = list())
#pairs(metrics[, c("complexity","spikiness","entropy","Background","Reads")],col=metrics$quality_manual)

################## Training models ##################

control <- trainControl(method="cv", number=10,classProbs = T)
metric <- "Accuracy"

library(DMwR)


set.seed(9560)
smote_train <- SMOTE(quality_manual ~ ., data  = trainData)
smote_outside <- train(quality_manual ~ ., data = smote_train,
					   method = "treebag",
					   nbagg = 50,
					   metric = "ROC",
					   trControl = control)
library(ROSE)
set.seed(9560)
rose_train <- ROSE::ROSE(quality_manual ~ ., data  = trainData)$data
set.seed(5627)
rose_outside <- train(quality_manual ~ ., data = rose_train,
					  method = "treebag",
					  nbagg = 50,
					  metric = "ROC",
					  trControl = control)
outside_models <- list(SMOTE = smote_outside,ROSE=rose_outside)
outside_resampling <- resamples(outside_models)


test_roc <- function(model, data) {
	library(pROC)
	roc_obj <- roc(data$quality_manual,
				   predict(model, data, type = "prob")[, "good"],
				   levels = c("good", "poor"))
	ci(roc_obj)
}

outside_test <- lapply(outside_models, test_roc, data = trainData)
outside_test <- lapply(outside_test, as.vector)
outside_test <- do.call("rbind", outside_test)
colnames(outside_test) <- c("lower", "ROC", "upper")
outside_test <- as.data.frame(outside_test)
summary(outside_resampling)


#Linear (lda), non-linear(CART,kNN), advanced
models=c("lda","rpart","knn","svmRadial","rf","gbm")
#model=models[1]
for (model in models){
	set.seed(7)
	fit <- train(quality_manual~., data=trainData, method=model, metric=metric, trControl=control)
	assign(paste0("fit.",model),fit)
	#results=resamples(list(get(paste0("fit.",model))))
}
#trainData=select(trainData,-c( Reads_per_Mb, Mode_GC, Mode_insert_size ,Sequencing_for_0.05x_cov))
#testData=select(testData,-c( Reads_per_Mb, Mode_GC, Mode_insert_size ,Sequencing_for_0.05x_cov))
#set.seed(7)
fit.gbm =  train(quality_manual~., data=trainData, method="gbm", metric=metric, trControl=control)
fit.rf <- train(quality_manual~., data=trainData, method=, metric=metric, trControl=control)
results <- resamples(list(lda=fit.lda,rpart=fit.rpart,knn=fit.knn,svmRadial=fit.svmRadial,rf=fit.rf,gbm=fit.gbm))
summary(results)
dotplot(results)
#default_knn_mod = train(quality_manual ~ .,data = trainData,method = "knn",trControl = trainControl(method = "cv", number = 5),preProcess = c("center", "scale"),tuneGrid = expand.grid(k = seq(1, 101, by = 2)))
#fit.gbm.imp <- varImp(fit.gbm)
#plot(fit.rf.imp, main="Variable Importance with Random Forest")
plot(fit.gbm, main="Model Accuracies with Random Forest")
plot(fit.rf, main="Model Accuracies with Random Forest")

#Running best model on testData/test set
predictions <- predict(fit.rf, na.omit(testData))
confusionMatrix(predictions, testData$quality_manual)

predictions <- predict(rose_outside, na.omit(testData))
confusionMatrix(predictions, testData$quality_manual)