args <- commandArgs(trailingOnly=TRUE)

CWD <- args[[1]] #project path
sample_list <- args[[2]] 
threads <- as.numeric(args[[3]])

dir.create(paste0(CWD,"Results/prediction_models"))

prediction_models_dir <- paste0(CWD,"Results/prediction_models/")


## ==============================
## Logging coloured (R)
## ==============================

use_color <- interactive() || Sys.getenv("TERM") != ""

if (use_color) {
  COL_INFO  <- "\033[34m"
  COL_OK    <- "\033[32m"
  COL_WARN  <- "\033[33m"
  COL_ERR   <- "\033[31m"
  COL_BOLD  <- "\033[1m"
  COL_RESET <- "\033[0m"
  
  SYM_OK   <- "✔"
  SYM_WARN <- "⚠"
  SYM_ERR  <- "✖"
} else {
  COL_INFO <- COL_OK <- COL_WARN <- COL_ERR <- COL_BOLD <- COL_RESET <- ""
  SYM_OK   <- "[OK]"
  SYM_WARN <- "[WARN]"
  SYM_ERR  <- "[ERROR]"
}

msg_info <- function(x) {
  cat(COL_INFO, COL_BOLD, x, COL_RESET, "\n", sep = "")
}

msg_ok <- function(x) {
  cat(COL_OK, COL_BOLD, SYM_OK, " ",  x, COL_RESET, "\n", sep = "")
}

msg_warn <- function(x) {
  cat(COL_WARN, COL_BOLD, SYM_WARN, " ",  x, COL_RESET, "\n", sep = "")
}

msg_error <- function(x) {
  cat(COL_ERR, COL_BOLD, SYM_ERR, " ",  x, COL_RESET, "\n", sep = "")
}


library(mlbench)
library(caret)
library(readxl)
library(caretEnsemble)
library(randomForest)
library(MLmetrics)
library(doParallel)

#Importing normalized matrix
data <- read.csv(paste0(CWD,"/Results/DEA_results/expression_matrix.csv"),sep=",",header=T)

#Importing metadata
samplesheet <- read.table(paste0(CWD,"/",sample_list), header = T, sep="\t")

rownames <- data$V1

data <- data[,-1]
rownames(data) <- rownames

#Assigning condition vector
if (!identical(colnames(data),samplesheet$sample)){
  msg_error("Sample-ids order in expression matriz and metadata does not match")
  stop("The condition is not met.", call. = FALSE)
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

msg_info("[1/2] Running RandomForest...")

set.seed(123)


rf_model <- train(condition ~ ., data=data, method="rf",metric="Accuracy",trControl=control)


msg_info("[2/2] Running GLMnet...")

set.seed(123)
glmnet_model <- train(condition ~ ., data=data, method="glmnet",metric="Accuracy",trControl=control)
}else{
  
msg_info("[1/2] Running RandomForest...")

set.seed(123)


rf_model <- train(condition ~ ., data=data, method="rf",metric="ROC",trControl=control)


msg_info("[2/2] Running GLMnet...")

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

#RANDOM FOREST

pred_rf <- rf_model$pred
pred$obs <- droplevels(pred_rf$obs)
clases_rf <- levels(pred_rf$obs)

if (length(clases_rf) == 2) {
  
  ## --------- BINARY CASE ---------
  positiva <- clases_rf[2]  # By default, the second class is the POSITIVE
  
  roc_rf <- roc(
    response  = pred_rf$obs,
    predictor = pred_rf[[positiva]],  
    levels    = clases_rf           
  )
  
  plot(roc_rf, main = "ROC curve - Random Forest")
  auc(roc_rf)
  
  pdf(paste0(prediction_models_dir,"/ROC-rf.pdf"), width = 10, height = 8)
  plot(roc_rf, main = "ROC curve - Random Forest")
  dev.off()
  
  
} else {
  
  ## --------- MULTICLASS CASE ---------
   # 1) Model global AUC
  prob_mat_rf <- as.matrix(pred_rf[, clases_rf])
  roc_multi_rf <- multiclass.roc(
    response  = pred_rf$obs,
    predictor = prob_mat_rf
  )
  print(roc_multi_rf$auc)   
  
  # 2) One vs all curves
  roc_list_rf <- lapply(clases_rf, function(cl) {
    resp_bin_rf <- factor(ifelse(pred_rf$obs == cl, cl, "other"),
                       levels = c("other", cl))  # cl es la positiva
    roc(
      response  = resp_bin_rf,
      predictor = pred[[cl]],
      levels    = c("other", cl)
    )
  })
  names(roc_list_rf) <- clases_rf
  
  # 3) Saving PDFs with all the curves
  pdf("ROC_RandomForest_multiclass.pdf")
  plot(
    roc_list[[1]],
    main = paste0("Multiclass ROC - Random Forest (AUC global = ",
                  round(as.numeric(roc_multi_rf$auc), 3), ")"),
    col  = 1
  )
  if (length(clases_rf) > 1) {
    for (i in 2:length(clases_rf)) {
      plot(roc_list_rf[[i]], add = TRUE, col = i)
    }
  }
  legend(
    "bottomright",
    legend = paste0(clases_rf, " (AUC=", round(sapply(roc_list_rf, auc), 3), ")"),
    col = seq_along(clases_rf),
    lty = 1,
    cex = 0.8
  )
  dev.off()  
  
  
}


#GLMNET

pred_glm <- glmnet_model$pred
pred_glm$obs <- droplevels(pred_glm$obs)
clases_glm <- levels(pred_glm$obs)


if (length(clases) == 2) {
  ## --------- BINARY CASE ---------
  positiva <- clases_glm[2]  
  
  roc_glmnet <- roc(
    response  = pred_glm$obs,
    predictor = pred_glm[[positiva]],  
    levels    = clases_glm             
  )
  
  plot(roc_glmnet, main = "ROC curve - GLMnet")
  auc(roc_glmnet)
  
  pdf(paste0(prediction_models_dir,"/ROC-glmnet.df"), width = 10, height = 8)
  plot(roc_glmnet, main = "ROC curve - GLMnet")
  dev.off()
  
} else {
  ## --------- MULTICLASS CASE ---------
  
   # 1) Model global AUC
  prob_mat_glm <- as.matrix(pred_glm[, clases_glm])
  roc_multi_glm <- multiclass.roc(
    response  = pred_glm$obs,
    predictor = prob_mat_glm
  )
  print(roc_multi_glm$auc)   
  
  # 2) One vs all curves
  roc_list_glm <- lapply(clases, function(cl) {
    resp_bin_glm <- factor(ifelse(pred_glm$obs == cl, cl, "other"),
                       levels = c("other", cl))
    roc(
      response  = resp_bin_glm,
      predictor = pred_glm[[cl]],
      levels    = c("other", cl)
    )
  })
  names(roc_list_glm) <- clases
  
  # 3) Saving PDFs for all the curves
  pdf("ROC_GLMnet_multiclass.pdf")
  plot(
    roc_list_glm[[1]],
    main = paste0("Multiclass ROC - GLMnet (AUC global = ",
                  round(as.numeric(roc_multi_glm$auc), 3), ")"),
    col  = 1
  )
  if (length(clases_glm) > 1) {
    for (i in 2:length(clases_glm)) {
      plot(roc_list_glm[[i]], add = TRUE, col = i)
    }
  }
  legend(
    "bottomright",
    legend = paste0(clases_glm, " (AUC=", round(sapply(roc_list_glm, auc), 3), ")"),
    col = seq_along(clases_glm),
    lty = 1,
    cex = 0.8
  )
  dev.off() 
}

