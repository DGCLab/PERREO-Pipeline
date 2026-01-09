args <- commandArgs(trailingOnly=TRUE)

CWD <- args[[1]] #project path
sample_list <- args[[2]] 
threads <- as.numeric(args[[3]])
positive_class <- args[[4]]


if (!dir.exists(paste0(CWD,"/prediction_models"))){
  dir.create(paste0(CWD,"/prediction_models"))}

prediction_models_dir <- paste0(CWD,"/prediction_models/")


suppressPackageStartupMessages({
  library(caret)
  library(pROC)
  library(dplyr)
  library(readr)
  library(tibble)
  library(doParallel)
})

dir.create(prediction_models_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(123)

#Importing normalized matrix
data <- read.csv(paste0(CWD,"/Results/DEA_results/expression_matrix.csv"),sep=",",header=T)

data <- data[!duplicated(data$X),]

rownames(data) <- data$X

data <- data[,-1]

#Importing metadata
samplesheet <- read.table(paste0(CWD,"/",sample_list), header = T, sep="\t")

# -------------------------
# 1) Filtering samples for sample_list
# -------------------------
samples_keep <- if ("sample" %in% colnames(samplesheet)) samplesheet$sample else samplesheet[[1]]
samples_keep <- samples_keep[!is.na(samples_keep) & samples_keep != ""]

samplesheet <- samplesheet %>%
  mutate(sample = as.character(sample),
         condition = as.factor(condition),
         strandedness = as.character(strandedness)) %>%
  filter(sample %in% samples_keep)

stopifnot(nrow(samplesheet) >= 4)

levels(samplesheet$condition) <- make.names(levels(samplesheet$condition))

# -------------------------
# 2) Building X (samples x features)
# -------------------------
samples <- samplesheet$sample
stopifnot(all(samples %in% colnames(data)))

X <- t(as.matrix(data[, samples, drop = FALSE]))
mode(X) <- "numeric"
X <- as.data.frame(X)
rownames(X) <- samples

y <- samplesheet$condition
names(y) <- samples

# removing features near-zero variance
nzv <- caret::nearZeroVar(X)
if (length(nzv) > 0) X <- X[, -nzv, drop = FALSE]

# -------------------------
# 3) Auto sampling
# -------------------------
pick_sampling_auto <- function(y) {
  tab <- sort(table(y), decreasing = TRUE)
  major_n <- as.numeric(tab[1])
  minor_n <- as.numeric(tab[length(tab)])
  ratio <- minor_n / major_n

  if (ratio >= 0.70) return("none")
  if (minor_n < 20) return("down")

  if (requireNamespace("DMwR", quietly = TRUE) || requireNamespace("DMwR2", quietly = TRUE)) return("smote")
  if (requireNamespace("ROSE", quietly = TRUE)) return("rose")
  return("down")
}
sampling_method <- pick_sampling_auto(y)
sampling_method <- if (sampling_method == "none") NULL else sampling_method

message("[INFO] sampling: ", ifelse(is.null(sampling_method), "none", sampling_method),
        " | dist: ", paste(names(table(y)), table(y), sep="=", collapse="; "))

# -------------------------
# 4) Split train/test
# -------------------------
idx <- caret::createDataPartition(y, p = 0.8, list = FALSE)
x_train <- X[idx, , drop = FALSE]
y_train <- y[idx]

x_test <- X[-idx, , drop = FALSE]
y_test <- y[-idx]

is_binary <- nlevels(y_train) == 2

if (is_binary) {
  pc <- make.names(positive_class)
  if (!pc %in% levels(y_train)) {
    stop("positive_class no coincide con niveles de condition: ", paste(levels(y_train), collapse=", "))
  }
  y_train <- factor(y_train, levels = c(pc, setdiff(levels(y_train), pc)))
  y_test  <- factor(y_test,  levels = levels(y_train))
} else {
  message("[INFO] Multiclass: se exporta ROC one-vs-rest por clase (positive_class no aplica).")
}

# -------------------------
# 5) Parallel
# -------------------------
threads <- max(1, threads)
cl <- parallel::makePSOCKcluster(threads)
doParallel::registerDoParallel(cl)
on.exit({ try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)

# -------------------------
# 6) trainControl
# -------------------------
if (is_binary) {
  ctrl <- trainControl(
    method = "repeatedcv", number = 5, repeats = 3,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final",
    sampling = sampling_method,
    allowParallel = TRUE
  )
  metric <- "ROC"
} else {
  ctrl <- trainControl(
    method = "repeatedcv", number = 5, repeats = 3,
    classProbs = TRUE,
    summaryFunction = multiClassSummary,
    savePredictions = "final",
    sampling = sampling_method,
    allowParallel = TRUE
  )
  metric <- "Accuracy"
}

# -------------------------
# 7) Training models
# -------------------------
glmnet_grid <- expand.grid(
  alpha = seq(0, 1, length.out = 6),
  lambda = 10^seq(-4, 1, length.out = 30)
)

glmnet_fit <- train(
  x = x_train, y = y_train,
  method = "glmnet",
  metric = metric,
  trControl = ctrl,
  tuneGrid = glmnet_grid,
  preProcess = c("center", "scale")
)

rf_grid <- expand.grid(mtry = unique(pmax(1, round(c(sqrt(ncol(x_train)), ncol(x_train)/3, ncol(x_train)/2)))))

rf_fit <- train(
  x = x_train, y = y_train,
  method = "rf",
  metric = metric,
  trControl = ctrl,
  tuneGrid = rf_grid
)

saveRDS(glmnet_fit, file.path(prediction_models_dir, "model_glmnet.rds"))
saveRDS(rf_fit,     file.path(prediction_models_dir, "model_rf.rds"))

write.csv(glmnet_fit$results, file.path(prediction_models_dir, "cv_results_glmnet.csv"), row.names = FALSE)
write.csv(rf_fit$results,     file.path(prediction_models_dir, "cv_results_rf.csv"), row.names = FALSE)

# -------------------------
# 8) Export predictions + ROC + metrics
# -------------------------
ovr_roc_curves <- function(truth, prob_df) {
  classes <- levels(truth)
  out <- list()
  for (cls in classes) {
    bin_truth <- factor(ifelse(truth == cls, cls, paste0("not_", cls)),
                        levels = c(cls, paste0("not_", cls)))
    roc_obj <- pROC::roc(bin_truth, prob_df[[cls]],
                         levels = rev(levels(bin_truth)), direction = "<", quiet = TRUE)
    out[[cls]] <- tibble(
      class = cls,
      fpr = 1 - roc_obj$specificities,
      tpr = roc_obj$sensitivities,
      threshold = roc_obj$thresholds,
      auc = as.numeric(pROC::auc(roc_obj))
    )
  }
  bind_rows(out)
}

eval_and_export <- function(fit, model_name) {
  prob_df <- predict(fit, newdata = x_test, type = "prob") %>% as.data.frame()
  pred    <- predict(fit, newdata = x_test, type = "raw")

  pred_out <- samplesheet %>%
    filter(sample %in% rownames(x_test)) %>%
    mutate(truth = as.character(y_test),
           pred  = as.character(pred)) %>%
    bind_cols(as_tibble(prob_df))

  write.csv(pred_out, file.path(prediction_models_dir, paste0("test_predictions_", model_name, ".csv")),
            row.names = FALSE)

  cm <- confusionMatrix(pred, y_test)

  base <- tibble(
    model = model_name,
    accuracy = unname(cm$overall["Accuracy"]),
    kappa = unname(cm$overall["Kappa"])
  )

  if (is_binary) {
    pos <- levels(y_test)[1]
    roc_obj <- pROC::roc(y_test, prob_df[[pos]], levels = rev(levels(y_test)),
                         direction = "<", quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))
    coords_best <- pROC::coords(roc_obj, x="best", best.method="youden",
                                ret=c("threshold","sensitivity","specificity"), transpose=FALSE)

    roc_df <- tibble(
      model = model_name, class = pos,
      fpr = 1 - roc_obj$specificities,
      tpr = roc_obj$sensitivities,
      threshold = roc_obj$thresholds,
      auc = auc_val
    )
    write.csv(roc_df, file.path(prediction_models_dir, paste0("roc_curve_", model_name, ".csv")),
              row.names = FALSE)

    p <- unname(cm$byClass["Pos Pred Value"])
    r <- unname(cm$byClass["Sensitivity"])
    f1 <- if ((p + r) > 0) 2*p*r/(p+r) else NA_real_

    bind_cols(base, tibble(
      auc = auc_val,
      sensitivity = unname(cm$byClass["Sensitivity"]),
      specificity = unname(cm$byClass["Specificity"]),
      precision = p,
      f1 = f1
    ))
  } else {
    roc_ovr <- ovr_roc_curves(y_test, prob_df) %>% mutate(model = model_name)
    write.csv(roc_ovr, file.path(prediction_models_dir, paste0("roc_ovr_", model_name, ".csv")),
              row.names = FALSE)

    mc_auc <- tryCatch({
      mroc <- pROC::multiclass.roc(response = y_test, predictor = as.matrix(prob_df), quiet = TRUE)
      as.numeric(mroc$auc)
    }, error = function(e) NA_real_)

    bind_cols(base, tibble(multiclass_auc = mc_auc))
  }
}

m1 <- eval_and_export(glmnet_fit, "glmnet")
m2 <- eval_and_export(rf_fit, "rf")

metrics_all <- bind_rows(m1, m2)
write.csv(metrics_all, file.path(prediction_models_dir, "metrics_summary.csv"), row.names = FALSE)

print(metrics_all)
cat("[OK] outputs en: ", normalizePath(prediction_models_dir), "\n", sep = "")

suppressPackageStartupMessages({
  library(pROC)
})

if (!dir.exists(prediction_models_dir)) {
  dir.create(prediction_models_dir, recursive = TRUE, showWarnings = FALSE)
}

.safe_roc <- function(truth_factor, probs_vec) {
  # Evita errores si hay NAs o si una clase no aparece en test
  ok <- is.finite(probs_vec) & !is.na(truth_factor)
  truth_factor <- droplevels(truth_factor[ok])
  probs_vec <- probs_vec[ok]
  if (nlevels(truth_factor) < 2) return(NULL)
  pROC::roc(truth_factor, probs_vec, levels = rev(levels(truth_factor)),
            direction = "<", quiet = TRUE)
}

if (exists("glmnet_fit") && exists("rf_fit") && exists("x_test") && exists("y_test") && exists("is_binary")) {

  if (isTRUE(is_binary)) {
    #Binary
    pos <- levels(y_test)[1]

    prob_glmnet <- tryCatch(predict(glmnet_fit, newdata = x_test, type = "prob")[[pos]],
                            error = function(e) NULL)
    prob_rf <- tryCatch(predict(rf_fit, newdata = x_test, type = "prob")[[pos]],
                        error = function(e) NULL)

    roc_glmnet <- if (!is.null(prob_glmnet)) .safe_roc(y_test, prob_glmnet) else NULL
    roc_rf     <- if (!is.null(prob_rf))     .safe_roc(y_test, prob_rf)     else NULL

    out_pdf <- file.path(prediction_models_dir, "ROC_test.pdf")
    pdf(out_pdf, width = 6, height = 6)

    if (!is.null(roc_glmnet)) {
      auc_glmnet <- as.numeric(pROC::auc(roc_glmnet))
      plot(roc_glmnet, main = "ROC (test)", legacy.axes = TRUE, col="blue")
    } else {
      plot.new()
      title("ROC (test)")
      text(0.5, 0.5, "No se pudo calcular ROC para GLMNET")
    }

    if (!is.null(roc_rf)) {
      auc_rf <- as.numeric(pROC::auc(roc_rf))
      plot(roc_rf, add = TRUE, col="red")
    } else {
      auc_rf <- NA_real_
    }

    if (!is.null(roc_glmnet) || !is.null(roc_rf)) {
      legend("bottomright",
             legend = c(
               if (!is.null(roc_glmnet)) paste0("GLMNET AUC=", round(auc_glmnet, 3)) else "GLMNET AUC=NA",
               if (!is.null(roc_rf))     paste0("RF AUC=", round(auc_rf, 3))         else "RF AUC=NA"
             ),
             lwd = 2, bty = "n")
    }
    dev.off()

    message("[OK] ROC PDF guardado: ", out_pdf)

  } else {
    # Multiclass: one-vs-rest
    .plot_ovr_pdf <- function(fit, model_name) {
      prob_df <- tryCatch(predict(fit, newdata = x_test, type = "prob") %>% as.data.frame(),
                          error = function(e) NULL)
      if (is.null(prob_df)) return(invisible(FALSE))

      classes <- levels(y_test)
      out_pdf <- file.path(prediction_models_dir, paste0("ROC_ovr_test_", model_name, ".pdf"))
      pdf(out_pdf, width = 7, height = 7)


      first <- TRUE
i <- 1  

for (cls in classes) {
  if (!cls %in% colnames(prob_df)) next

  bin_truth <- factor(
    ifelse(y_test == cls, cls, paste0("not_", cls)),
    levels = c(cls, paste0("not_", cls))
  )

  roc_obj <- .safe_roc(bin_truth, prob_df[[cls]])
  if (is.null(roc_obj)) next

  col_i <- base_colors[i]

  if (first) {
    plot(
      roc_obj,
      main = paste0("ROC one-vs-rest (test) - ", model_name),
      legacy.axes = TRUE,
      col = col_i,
      lwd = 2
    )
    first <- FALSE
  } else {
    plot(roc_obj, add = TRUE, col = col_i, lwd = 2)
  }

  i <- i + 1
}

      if (first) {
        plot.new()
        title(paste0("ROC one-vs-rest (test) - ", model_name))
        text(0.5, 0.5, "No se pudieron calcular curvas ROC (clases insuficientes o datos)")
      } else {
    legend("bottomright",
       legend = classes,
       col = base_colors[seq_along(classes)],
       lwd = 2)

      }

      dev.off()
      message("[OK] ROC OVR PDF guardado: ", out_pdf)
      invisible(TRUE)
    }

    .plot_ovr_pdf(glmnet_fit, "glmnet")
    .plot_ovr_pdf(rf_fit, "rf")
  }

} else {
  message("[WARN] No se generaron PDFs ROC: faltan objetos requeridos (glmnet_fit/rf_fit/x_test/y_test/is_binary).")
}

# ---- Export top variables (RF + GLMNET) ----
top_n <- 50  

get_top <- function(fit, n = 50) {
  imp <- caret::varImp(fit)$importance
  score <- if ("Overall" %in% colnames(imp)) imp$Overall else rowMeans(imp[, sapply(imp, is.numeric), drop = FALSE], na.rm = TRUE)
  out <- data.frame(feature = rownames(imp), importance = as.numeric(score))
  out <- out[order(out$importance, decreasing = TRUE), ]
  head(out, n)
}

write.csv(get_top(rf_fit, top_n),
          file.path(prediction_models_dir, paste0("top", top_n, "_features_rf.csv")),
          row.names = FALSE)

write.csv(get_top(glmnet_fit, top_n),
          file.path(prediction_models_dir, paste0("top", top_n, "_features_glmnet.csv")),
          row.names = FALSE)


