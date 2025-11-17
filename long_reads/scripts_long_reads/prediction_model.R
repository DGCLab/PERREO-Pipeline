args <- commandArgs(trailingOnly=TRUE)

CWD <- args[[1]] #project path
sample_list <- args[[2]] 

dir.create(paste0(CWD,"/prediction_models"))

prediction_models_dir <- paste0(CWD,"/prediction_models/")


library(mlbench)
library(caret)
library(readxl)
library(caretEnsemble)
library(randomForest)
library(MLmetrics)
#Importing normalized matrix
data <- read.csv(paste0(CWD,"/SAMPLES/DEA_results/expression_matrix.csv"),sep=",",header=T)

#Importing metadata
samplesheet <- read.table(paste0(CWD,"/",sample_list), header = T, sep="\t")

rownames <- data$V1

data <- data[,-1]
rownames(data) <- rownames

#Assigning condition vector
if (!identical(colnames(data),samplesheet$sample)){
  stop("Sample-ids order in expression matriz and metadata does not match")
}

data <- t(data)
condition <- samplesheet$condition
data <- cbind(data,condition)
data <- as.data.frame(data)
#data$condition <- factor(data$condition, levels = levels(data$condition))

if (length(unique(condition))){

control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats=3,
  classProbs = TRUE,
  savePredictions = "final",
  summaryFunction = multiClassSummary)
}else{
control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats=3,
  classProbs = TRUE,
  savePredictions = "final",
  summaryFunction = twoClassSummary)
}



#models <- c("rf","svmRadial","gbm","glmnet","lda")

#library(doParallel)
#registerDoSEQ()  # Deactivating paralelism


set.seed(123)
rf_model <- train(condition ~ ., data=data, method="rf",trControl=control)

set.seed(123)
svmradial_model <- train(condition ~ ., data=data, method="svmRadial",trControl=control)

set.seed(123)
gbm_model <- train(condition ~ ., data=data, method="gbm",trControl=control)

set.seed(123)
glmnet_model <- train(condition ~ ., data=data, method="glmnet",trControl=control)


#Comparing results of the 4 models

results <- resamples(list(RF=rf_model,SVMR=svmradial_model,GBM=gbm_model,
                          GLMnet=glmnet_model))

#Generating plots 

summary(results)
p1<-bwplot(results)
p2<-dotplot(results)
p3<-densityplot(results)


#Saving plots 
png(paste0(prediction_models_dir,"bwplot.png"), width = 1200, height = 800, res = 150)
print(p1)
dev.off()
png(paste0(prediction_models_dir,"dotplot.png"), width = 1200, height = 800, res = 150)
print(p2)
dev.off()
png(paste0(prediction_models_dir,"density.png"), width = 1200, height = 800, res = 150)
print(p3)
dev.off()


#Extracting metrics in table format

results_table <- results$values

#Calculating means for model 

means <- data.frame(
  Model = gsub("~ROC|~Sens|~Spec", "", names(results_table)[-1]),
  Metric = gsub(".*~", "", names(results_table)[-1]),
  Mean = colMeans(results_table[,-1])
)

write.csv(means,paste0(prediction_models_dir,"models_metrics.csv"), row.names = FALSE)

#Cleaner format

library(tidyr)
results_wide <- means %>%
  pivot_wider(names_from = Metric, values_from = Mean)

write.csv(results_wide,paste0(prediction_models_dir,"summary_metrics.csv"), row.names = FALSE)

#ROC curves

library(pROC)

head(rf_model$pred)
roc_rf <- roc(rf_model$pred$obs, rf_model$pred$LUAD)
r_rf<-plot(roc_rf, main = "ROC curve - Random Forest")
auc(roc_rf)

head(svmradial_model$pred)
roc_svmr <- roc(svmradial_model$pred$obs, svmradial_model$pred$LUAD)
r_svmr<-plot(roc_svmr, main = "ROC curve - SVMRadial")
auc(roc_svmr)

head(gbm_model$pred)
roc_gbm <- roc(gbm_model$pred$obs, gbm_model$pred$LUAD)
r_gbm<-plot(roc_gbm, main = "ROC curve - GBM")
auc(roc_gbm)

head(glmnet_model$pred)
roc_glmnet <- roc(glmnet_model$pred$obs, glmnet_model$pred$LUAD)
r_glmnet<-plot(roc_glmnet, main = "ROC curve - GLMnet")
auc(roc_glmnet)

#Saving ROC curves
#Saving plots 
png(paste0(prediction_models_dir,"ROC-rf.png"), width = 1200, height = 800, res = 150)
plot(roc_rf, main = "ROC curve - Random Forest")
dev.off()
png(paste0(prediction_models_dir,"ROC-svmr.png"), width = 1200, height = 800, res = 150)
plot(roc_svmr, main = "ROC curve - SVMRadial")
dev.off()
png(paste0(prediction_models_dir,"ROC-gbm.png"), width = 1200, height = 800, res = 150)
plot(roc_gbm, main = "ROC curve - GBM")
dev.off()
png(paste0(prediction_models_dir,"ROC-glmnet.png"), width = 1200, height = 800, res = 150)
plot(roc_glmnet, main = "ROC curve - GLMnet")
dev.off()




