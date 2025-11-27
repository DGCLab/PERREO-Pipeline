args <- commandArgs(trailingOnly=TRUE)

CWD <- args[[1]] #project path
sample_list <- args[[2]] 
threads <- as.numeric(args[[3]])

dir.create(paste0(CWD,"/prediction_models"))

prediction_models_dir <- paste0(CWD,"/prediction_models/")


library(mlbench)
library(caret)
library(readxl)
library(caretEnsemble)
library(randomForest)
library(MLmetrics)
library(doParallel)

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


num_cores <- threads - 1
cl <- makePSOCKcluster(num_cores)

registerDoParallel(cl)

if (length(unique(condition))>2){

control <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats=5,
  classProbs = TRUE,
  savePredictions = "final",
  summaryFunction = multiClassSummary,
allowParallel=TRUE)
}else{
control <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats=5,
  classProbs = TRUE,
  savePredictions = "final",
  summaryFunction = twoClassSummary,
allowParallel=TRUE)
}


if (length(unique(condition))>2){

print("Running RandomForest...")

set.seed(123)


rf_model <- train(condition ~ ., data=data, method="rf",metric="Accuracy",trControl=control)


print("Running GLMnet...")

set.seed(123)
glmnet_model <- train(condition ~ ., data=data, method="glmnet",metric="Accuracy",trControl=control)
}else{
print("Running RandomForest...")

set.seed(123)


rf_model <- train(condition ~ ., data=data, method="rf",metric="ROC",trControl=control)


print("Running GLMnet...")

set.seed(123)
glmnet_model <- train(condition ~ ., data=data, method="glmnet",metric="ROC",trControl=control)
}


#Comparing results of the 4 models

results <- resamples(list(RF=rf_model,
                          GLMnet=glmnet_model))

stopCluster(cl)
registerDoSEQ()

#Generating plots 

summary(results)
p1<-bwplot(results)
p2<-dotplot(results)
p3<-densityplot(results)


#Saving plots 
pdf(paste0(prediction_models_dir,"/bwplot.pdf"), width = 10, height = 8)
print(p1)
dev.off()
pdf(paste0(prediction_models_dir,"/dotplot.pdf"), width = 10, height = 8)
print(p2)
dev.off()
pdf(paste0(prediction_models_dir,"/density.pdf"), width = 10, height = 8)
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

write.csv(means,paste0(prediction_models_dir,"/models_metrics.csv"), row.names = FALSE)

#Cleaner format

library(tidyr)
results_wide <- means %>%
  pivot_wider(names_from = Metric, values_from = Mean)

write.csv(results_wide,paste0(prediction_models_dir,"/summary_metrics.csv"), row.names = FALSE)

#ROC curves

library(pROC)

pred <- rf_model$pred
pred$obs <- droplevels(pred$obs)
clases <- levels(pred$obs)


if (length(clases) == 2) {
  
  ## --------- BINARY CASE ---------
  positiva <- clases[2]  # By default, the second class is the POSITIVE
  
  roc_rf <- roc(
    response  = pred$obs,
    predictor = pred[[positiva]],  
    levels    = clases           
  )
  
  plot(roc_rf, main = "ROC curve - Random Forest")
  auc(roc_rf)
  
  pdf(paste0(prediction_models_dir,"/ROC-rf.pdf"), width = 10, height = 8)
  plot(roc_rf, main = "ROC curve - Random Forest")
  dev.off()
  
  
} else {
  
  ## --------- MULTICLASS CASE ---------
   # 1) Model global AUC
  prob_mat <- as.matrix(pred[, clases])
  roc_multi <- multiclass.roc(
    response  = pred$obs,
    predictor = prob_mat
  )
  print(roc_multi$auc)   
  
  # 2) One vs all curves
  roc_list <- lapply(clases, function(cl) {
    resp_bin <- factor(ifelse(pred$obs == cl, cl, "other"),
                       levels = c("other", cl))  # cl es la positiva
    roc(
      response  = resp_bin,
      predictor = pred[[cl]],
      levels    = c("other", cl)
    )
  })
  names(roc_list) <- clases
  
  # 3) Saving PDFs with all the curves
  pdf("ROC_RandomForest_multiclass.pdf")
  plot(
    roc_list[[1]],
    main = paste0("Multiclass ROC - Random Forest (AUC global = ",
                  round(as.numeric(roc_multi$auc), 3), ")"),
    col  = 1
  )
  if (length(clases) > 1) {
    for (i in 2:length(clases)) {
      plot(roc_list[[i]], add = TRUE, col = i)
    }
  }
  legend(
    "bottomright",
    legend = paste0(clases, " (AUC=", round(sapply(roc_list, auc), 3), ")"),
    col = seq_along(clases),
    lty = 1,
    cex = 0.8
  )
  dev.off()  
  
  
}


pred <- glmnet_model$pred
pred$obs <- droplevels(pred$obs)
clases <- levels(pred$obs)


if (length(clases) == 2) {
  ## --------- BINARY CASE ---------
  positiva <- clases[2]  
  
  roc_glmnet <- roc(
    response  = pred$obs,
    predictor = pred[[positiva]],  
    levels    = clases             
  )
  
  plot(roc_glmnet, main = "ROC curve - Random Forest")
  auc(roc_glmnet)
  
  pdf(paste0(prediction_models_dir,"/ROC-rf.df"), width = 10, height = 8)
  plot(roc_glmnet, main = "ROC curve - Random Forest")
  dev.off()
  
} else {
  ## --------- MULTICLASS CASE ---------
  
   # 1) Model global AUC
  prob_mat <- as.matrix(pred[, clases])
  roc_multi <- multiclass.roc(
    response  = pred$obs,
    predictor = prob_mat
  )
  print(roc_multi$auc)   
  
  # 2) One vs all curves
  roc_list <- lapply(clases, function(cl) {
    resp_bin <- factor(ifelse(pred$obs == cl, cl, "other"),
                       levels = c("other", cl))
    roc(
      response  = resp_bin,
      predictor = pred[[cl]],
      levels    = c("other", cl)
    )
  })
  names(roc_list) <- clases
  
  # 3) Saving PDFs for all the curves
  pdf("ROC_RandomForest_multiclass.pdf")
  plot(
    roc_list[[1]],
    main = paste0("Multiclass ROC - Random Forest (AUC global = ",
                  round(as.numeric(roc_multi$auc), 3), ")"),
    col  = 1
  )
  if (length(clases) > 1) {
    for (i in 2:length(clases)) {
      plot(roc_list[[i]], add = TRUE, col = i)
    }
  }
  legend(
    "bottomright",
    legend = paste0(clases, " (AUC=", round(sapply(roc_list, auc), 3), ")"),
    col = seq_along(clases),
    lty = 1,
    cex = 0.8
  )
  dev.off() 
}
