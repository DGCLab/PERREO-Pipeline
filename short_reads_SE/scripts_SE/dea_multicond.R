###############################################################################
####                                                                       ####
####                      DIFFERENTIAL EXPRESSION ANALYSIS                 ####
####                                                                       ####
####             Authors:                                                  ####
####                                                                       ####
####                                                                       ####
###############################################################################


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


version <- "v1.0"
cat(sprintf("
           /^ ^\\                                    /^ ^\\
          / 0 0 \\       PERREO PIPELINE %s       / 0 0 \\
          V\\ Y /V  ------------------------------  V\\ Y /V
           / - \\     Processing Differential        / - \\
           |    \\      Expression Analysis         /    |
           ||(__V)                                (__V)||

------------------------------------------------------------

", version))



options(warn = -1)

args <- commandArgs(trailingOnly=TRUE)

batch <- args[[1]]
sample_list <- args[[2]]
method <- args[[3]]
CWD <- args[[4]]
repeatmasker_annotation_gtf <- args[[5]]
k_num <- as.numeric(args[[6]])
FDR_thr <- as.numeric(args[[7]])
log2FC_thr <- as.numeric(args[[8]])
SAMPLES_DIR <- paste0(CWD,"/Results/")
DEA_results_DIR <- paste0(SAMPLES_DIR,"DEA_results")


if (batch == "yes"){batch = TRUE}else{batch=FALSE}

## Load libraries
suppressPackageStartupMessages({
library(dplyr)
library(DESeq2)
library(readr)
library(ggplot2)
library(pheatmap)
library(ggpubr)
library(ggrepel)
library(grid)
library(rtracklayer)
library(edgeR)
library(limma)
library(EDASeq)
library(RUVSeq)
library(pdftools)
})
  
# Loading the metadata:

sample_list <- read.table(paste0(CWD,"/",sample_list), header = T, sep="\t")

## Prepare the data for DESeq2:

cts_raw <- read.table(paste0(SAMPLES_DIR,"count_data.txt"))
rownames <- rownames(cts_raw)

cts <- as.matrix(cts_raw)
rownames(cts) <- rownames
samples <- colnames(cts)

if (!identical(samples, sample_list$sample)) {
  sample_list <- sample_list[order(sample_list$sample),]
}

if (identical(as.character(samples), as.character(sample_list$sample))) {
  msg_ok("Order of samples ids in `cts` and `samplesheet` do match!")
} else {  
  msg_error("Order of samples ids in `cts` and `samplesheet` do not match!")
  stop("The condition is not met.", call. = FALSE)
}

condition <- sample_list$condition
coldata <- within(data.frame(row.names = samples,
                             condition = factor(condition)), {
                               if ("batch" %in% colnames(sample_list))
                                 batch <- sample_list$batch})

if (ncol(cts) != nrow(coldata)) {
  msg_error("Number of samples in `cts` and `coldata` do not match!")
  stop("The condition is not met.", call. = FALSE)
}

cts <- round(cts)

dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = ~  condition 
)

## Filter genes which have at least a zero between all samples

keep <- apply(counts(dds), 1, function(x) all(x > 0))
dds <- dds[keep, ]

mat <- counts(dds)

## Set all possible contrasts

colData(dds)$condition <- factor(colData(dds)$condition)
levels(colData(dds)$condition)

lv <- levels(colData(dds)$condition)
pairs <- combn(lv, 2, simplify = FALSE)

outdir <- paste0(DEA_results_DIR,"/Contrasts")
dir.create(outdir, showWarnings = FALSE)

## Correct Batch Effect

if(batch == TRUE){
  ruvg <- TRUE
  msg_info("Removing Batch Effect...")
  
  if ("batch" %in% colnames(sample_list)) {
    ruvg <- FALSE
  }
  
  if(ruvg == TRUE){
  dds0 <- DESeq(dds, test = "LRT", reduced = ~ 1)
  lrt  <- results(dds0)
  
  keep <- lrt$baseMean > 10 & !is.na(lrt$pvalue)
  pvec <- lrt$pvalue[keep]
  
  thr  <- quantile(pvec, 0.75, na.rm = TRUE)
  neg_controls <- rownames(lrt)[keep & lrt$pvalue >= thr]
  
  set <- newSeqExpressionSet(
    counts = mat,
    normalizedCounts = cpm(mat), # For PCA pre-BATCH
    phenoData = data.frame(coldata, row.names = colnames(cts)))
  
  
  res <- RUVg(x=set, cIdx=neg_controls, k=k_num ,isLog = F)
  mat.cpm.ruv <- res@assayData$normalizedCounts
  
  Wdf <- as.data.frame(pData(res)[, grep("^W_", colnames(pData(res))), drop = FALSE]) 
  
  # Add to metadata
  sample_list <- cbind(sample_list, Wdf)
  
  if (method == "edgeR") {

    msg_info("DEA is being performed by EdgeR")
    
    dge <- DGEList(counts = mat, samples = samples)
    dge <- calcNormFactors(dge, method="TMM")
    
    form <- as.formula(paste("~ ", paste(colnames(Wdf), collapse = " + "), "+ condition"))
    design <- model.matrix(form, data = sample_list)
    
    dge <- estimateDisp(dge, design)
    fit <- glmQLFit(dge, design)
    
    mat.dge <- cpm(dge) # For model prediction
    
    coef_names <- colnames(design)
    cond_coef <- setNames(rep(NA_character_, length(lv)), lv)
    
    for (lev in lv[-1]) {
      nm <- paste0("condition", lev)
      if (nm %in% coef_names) cond_coef[lev] <- nm
    }
    
    cond_coef[lv[1]] <- NA
    
    for (p in pairs) {
      den <- p[1]   
      num <- p[2]   
      cname <- paste0(num, "_vs_", den)
      message("Processing contrast: ", cname)
      
      # Construir contraste
      # Casos:
      # 1) num != baseline y den = baseline  -> coef = + condition<num>
      # 2) num = baseline y den != baseline  -> coef = - condition<den>
      # 3) num != baseline y den != baseline -> contrast vector: +cond<num> -cond<den>
      if (is.na(cond_coef[num]) && is.na(cond_coef[den])) {
        msg_error("Both conditions are baselines.")
        stop("The condition is not met.", call. = FALSE)
      }
      
      if (!is.na(cond_coef[num]) && is.na(cond_coef[den])) {
        # num vs baseline
        qlf <- glmQLFTest(fit, coef = which(coef_names == cond_coef[num]))
      } else if (is.na(cond_coef[num]) && !is.na(cond_coef[den])) {
        # baseline vs den  => -(den vs baseline)
        v <- rep(0, ncol(design))
        v[coef_names == cond_coef[den]] <- -1
        qlf <- glmQLFTest(fit, contrast = v)
      } else {
        # num vs den (ambos no-baseline)
        v <- rep(0, ncol(design))
        v[coef_names == cond_coef[num]] <-  1
        v[coef_names == cond_coef[den]] <- -1
        qlf <- glmQLFTest(fit, contrast = v)
      }
      
      res_df <- topTags(qlf, n = Inf)$table
      res_df$gene <- rownames(res_df)
      
      res_filtered <- res_df[abs(res_df$logFC) > log2FC_thr & res_df$FDR < FDR_thr, , drop = FALSE]
      assign(paste0("res_filtered","_contrast_", cname), res_filtered)
      
      write.csv(res_df, file = file.path(outdir, paste0("contrast_", cname, ".csv")), row.names = F)
      msg_ok(paste0("Saved: ", file.path(outdir, paste0("contrast_", cname, ".csv"))))
      
      assign(paste0("res","_contrast_", cname), res_df)
    }
    
  } else if (method == "DESeq2") {
    
    msg_info("DEA is being performed by DESeq2")
    
    form <- as.formula(paste("~ ", paste(colnames(Wdf), collapse = " + "), "+ condition"))
    
    dds <- DESeqDataSetFromMatrix(
      countData = mat,
      colData = sample_list,
      design = form)
      
    dds <- DESeq(dds)
    
    for (p in pairs) {
      num <- p[2]   
      den <- p[1]   
      cname <- paste0(num, "_vs_", den)
      message("Processing contrast: ", cname)
      
      
      res <- results(dds, contrast = c("condition", num, den))
      res_df <- as.data.frame(res)
      res_df$gene <- rownames(res_df)
      
      res_filtered <- subset(res_df, abs(res_df$log2FoldChange) > log2FC_thr & res_df$padj < FDR_thr) 
      assign(paste0("res_filtered","_contrast_", cname), res_filtered)
      
      write.csv(res_df, file = file.path(outdir, paste0("contrast_", cname, ".csv")), row.names = F)
      msg_ok(paste0("Saved: ", file.path(outdir, paste0("contrast_", cname, ".csv"))))
      
      assign(paste0("res","_contrast_", cname), res_df)
      
    }
    
    vsd <- varianceStabilizingTransformation(dds) # For model prediction
  }
  
  } else {
   if (method == "edgeR") {

    msg_info("DEA is being performed by EdgeR")
    
    dge <- DGEList(counts = mat, samples = samples)
    dge <- calcNormFactors(dge, method="TMM")
    
    design <- model.matrix(~ batch + condition, data = sample_list)
    dge <- estimateDisp(dge, design)
    fit <- glmQLFit(dge, design)
    
    mat.dge <- cpm(dge) # For model prediction
    
    coef_names <- colnames(design)
    cond_coef <- setNames(rep(NA_character_, length(lv)), lv)
    
    for (lev in lv[-1]) {
      nm <- paste0("condition", lev)
      if (nm %in% coef_names) cond_coef[lev] <- nm
    }
    
    cond_coef[lv[1]] <- NA
    
    for (p in pairs) {
      den <- p[1]   
      num <- p[2]   
      cname <- paste0(num, "_vs_", den)
      message("Processing contrast: ", cname)
      
      # Construir contraste
      # Casos:
      # 1) num != baseline y den = baseline  -> coef = + condition<num>
      # 2) num = baseline y den != baseline  -> coef = - condition<den>
      # 3) num != baseline y den != baseline -> contrast vector: +cond<num> -cond<den>
      if (is.na(cond_coef[num]) && is.na(cond_coef[den])) {
        msg_error("Both conditions are baselines.")
        stop("The condition is not met.", call. = FALSE)
      }
      
      if (!is.na(cond_coef[num]) && is.na(cond_coef[den])) {
        # num vs baseline
        qlf <- glmQLFTest(fit, coef = which(coef_names == cond_coef[num]))
      } else if (is.na(cond_coef[num]) && !is.na(cond_coef[den])) {
        # baseline vs den  => -(den vs baseline)
        v <- rep(0, ncol(design))
        v[coef_names == cond_coef[den]] <- -1
        qlf <- glmQLFTest(fit, contrast = v)
      } else {
        # num vs den (ambos no-baseline)
        v <- rep(0, ncol(design))
        v[coef_names == cond_coef[num]] <-  1
        v[coef_names == cond_coef[den]] <- -1
        qlf <- glmQLFTest(fit, contrast = v)
      }
      
      res_df <- topTags(qlf, n = Inf)$table
      res_df$gene <- rownames(res_df)
      
      res_filtered <- res_df[abs(res_df$logFC) > log2FC_thr & res_df$FDR < FDR_thr, , drop = FALSE]
      assign(paste0("res_filtered","_contrast_", cname), res_filtered)
      
      write.csv(res_df, file = file.path(outdir, paste0("contrast_", cname, ".csv")), row.names = F)
      msg_ok(paste0("Saved: ", file.path(outdir, paste0("contrast_", cname, ".csv"))))
      
      assign(paste0("res","_contrast_", cname), res_df)
      
    }
    
  } else if (method == "DESeq2") 
    {
    
    msg_info("DEA is being performed by DESeq2")
    
    dds <- DESeqDataSetFromMatrix(
      countData = mat,
      colData = sample_list,
      design = ~ batch + condition)
    
    dds <- DESeq(dds)
    
    for (p in pairs) {
      num <- p[2]   
      den <- p[1]   
      cname <- paste0(num, "_vs_", den)
      message("Processing contrast: ", cname)
      
      
      res <- results(dds, contrast = c("condition", num, den))
      res_df <- as.data.frame(res)
      res_df$gene <- rownames(res_df)
      
      res_filtered <- subset(res_df, abs(res_df$log2FoldChange) > log2FC_thr & res_df$padj < FDR_thr) 
      assign(paste0("res_filtered","_contrast_", cname), res_filtered)
      
      write.csv(res_df, file = file.path(outdir, paste0("contrast_", cname, ".csv")), row.names = F)
      msg_ok(paste0("Saved: ", file.path(outdir, paste0("contrast_", cname, ".csv"))))
      
      assign(paste0("res","_contrast_", cname), res_df)
      
    }
    
    vsd <- varianceStabilizingTransformation(dds) # For model prediction
    }
  }
  
  # Once the batch effect is corrected we use the matrix for visualization
  # PCA
  pca <- prcomp(t(set@assayData$normalizedCounts),scale. = T)
  pca_df <- as.data.frame(pca$x)
  pca_df$condition <- condition
  
  base_colors <- c("#804A45", "#455F80","#F5A553", "#BAC6D4", "#2E8B57", "#323840", "#FFF1C2")
  palette_auto <- base_colors[seq_len(length(unique(pca_df$condition)))]
  
  p <- ggplot(pca_df, aes(x=PC1, y=PC2, color=condition)) +
    geom_point(size=3) +
    scale_color_manual(values = palette_auto) +
    stat_ellipse(aes(group = condition), type = "norm", level = 0.95, linetype = "dashed") +
    labs(title="PCA - Before Batch Effect Corrected",
         x=paste0("PC1 (", round(100*summary(pca)$importance[2,1], 1), "% var)"),
         y=paste0("PC2 (", round(100*summary(pca)$importance[2,2], 1), "% var)")) +
    theme_minimal()
  
  ggsave(paste0(DEA_results_DIR,"/pca_nobatch.png"), plot = p, width = 8, height = 6, dpi = 300)
  
  
  pca_corrected <- prcomp(t(mat.cpm.ruv),scale. = T)
  pca_df_corrected <- as.data.frame(pca_corrected$x)
  pca_df_corrected$condition <- condition
  
  p_corrected <- ggplot(pca_df_corrected, aes(x=PC1, y=PC2, color=condition)) +
    geom_point(size=3) +
    scale_color_manual(values = palette_auto) +
    stat_ellipse(aes(group = condition), type = "norm", level = 0.95, linetype = "dashed") +
    labs(title="PCA - Batch Effect Corrected",
         x=paste0("PC1 (", round(100*summary(pca)$importance[2,1], 1), "% var)"),
         y=paste0("PC2 (", round(100*summary(pca)$importance[2,2], 1), "% var)")) +
    theme_minimal()
  
  ggsave(paste0(DEA_results_DIR,"/pca_corrected.png"), plot = p_corrected, width = 8, height = 6, dpi = 300)
  ggsave(paste0(DEA_results_DIR,"/pca_corrected.pdf"), plot = p_corrected, width = 8, height = 6, dpi = 300)

  } else {
  if (method == "edgeR") {
    msg_info("DEA is being performed by EdgeR")
    
    dge <- DGEList(counts = mat, samples = samples)
    dge <- calcNormFactors(dge, method="TMM")
    
    design <- model.matrix(~ condition, data = sample_list)
    dge <- estimateDisp(dge, design)
    fit <- glmQLFit(dge, design)
    
    mat.dge <- cpm(dge) # For model prediction & PCA
    
    coef_names <- colnames(design)
    cond_coef <- setNames(rep(NA_character_, length(lv)), lv)
    
    for (lev in lv[-1]) {
      nm <- paste0("condition", lev)
      if (nm %in% coef_names) cond_coef[lev] <- nm
    }
    
    cond_coef[lv[1]] <- NA
    
    for (p in pairs) {
      den <- p[1]   
      num <- p[2]   
      cname <- paste0(num, "_vs_", den)
      message("Processing contrast: ", cname)
      
      # Construir contraste
      # Casos:
      # 1) num != baseline y den = baseline  -> coef = + condition<num>
      # 2) num = baseline y den != baseline  -> coef = - condition<den>
      # 3) num != baseline y den != baseline -> contrast vector: +cond<num> -cond<den>
      if (is.na(cond_coef[num]) && is.na(cond_coef[den])) {
        msg_error("Error: both conditions are baselines.")
        stop("The condition is not met.", call. = FALSE)
      }
      
      if (!is.na(cond_coef[num]) && is.na(cond_coef[den])) {
        # num vs baseline
        qlf <- glmQLFTest(fit, coef = which(coef_names == cond_coef[num]))
      } else if (is.na(cond_coef[num]) && !is.na(cond_coef[den])) {
        # baseline vs den  => -(den vs baseline)
        v <- rep(0, ncol(design))
        v[coef_names == cond_coef[den]] <- -1
        qlf <- glmQLFTest(fit, contrast = v)
      } else {
        # num vs den (ambos no-baseline)
        v <- rep(0, ncol(design))
        v[coef_names == cond_coef[num]] <-  1
        v[coef_names == cond_coef[den]] <- -1
        qlf <- glmQLFTest(fit, contrast = v)
      }
      
      res_df <- topTags(qlf, n = Inf)$table
      res_df$gene <- rownames(res_df)
      
      
      res_filtered <- res_df[abs(res_df$logFC) > log2FC_thr & res_df$FDR < FDR_thr, , drop = FALSE]
      assign(paste0("res_filtered","_contrast_", cname), res_filtered)
      
      write.csv(res_df, file = file.path(outdir, paste0("contrast_", cname, ".csv")), row.names = F)
      msg_ok(paste0("Saved: ", file.path(outdir, paste0("contrast_", cname, ".csv"))))
      
      assign(paste0("results","_contrast_", cname), res_df)
    }
    
    pca <- prcomp(t(mat.dge),scale. = T)
    pca_df <- as.data.frame(pca$x)
    pca_df$condition <- condition
    
    base_colors <- c("#804A45", "#455F80","#F5A553", "#BAC6D4", "#2E8B57", "#323840", "#FFF1C2")
    palette_auto <- base_colors[seq_len(length(unique(pca_df$condition)))]
    
    
    p <- ggplot(pca_df, aes(x=PC1, y=PC2, color=condition)) +
      scale_color_manual(values = palette_auto) +
      geom_point(size=3) +
      stat_ellipse(aes(group = condition), type = "norm", level = 0.95, linetype = "dashed") +
      labs(title="PCA",
           x=paste0("PC1 (", round(100*summary(pca)$importance[2,1], 1), "% var)"),
           y=paste0("PC2 (", round(100*summary(pca)$importance[2,2], 1), "% var)")) +
      theme_minimal() +
      geom_point(size=3) +
      stat_ellipse(aes(group = condition), type = "norm", level = 0.95, linetype = "dashed") +
      labs(title="PCA",
           x=paste0("PC1 (", round(100*summary(pca)$importance[2,1], 1), "% var)"),
           y=paste0("PC2 (", round(100*summary(pca)$importance[2,2], 1), "% var)")) +
      theme_minimal()
    
    ggsave(paste0(DEA_results_DIR,"/pca.png"), plot = p, width = 8, height = 6, dpi = 300)
    ggsave(paste0(DEA_results_DIR,"/pca.pdf"), plot = p, width = 8, height = 6, dpi = 300)
    
    
  } else if (method == "DESeq2") 
  {
    
    msg_info("DEA is being performed by DESeq2")
    
    dds <- DESeqDataSetFromMatrix(
      countData = mat,
      colData = sample_list,
      design = ~ condition)
    
    dds <- DESeq(dds)
    
    vsd <- varianceStabilizingTransformation(dds) # For model prediction & PCA
    
    for (p in pairs) {
      num <- p[2]   
      den <- p[1]   
      cname <- paste0(num, "_vs_", den)
      message("Processing contrast: ", cname)
      
      
      res <- results(dds, contrast = c("condition", num, den))
      res_df <- as.data.frame(res)
      res_df$gene <- rownames(res_df)
      
      res_filtered <- subset(res_df, abs(log2FoldChange) > log2FC_thr & padj < FDR_thr) 
      assign(paste0("res_filtered","_contrast_", cname), res_filtered)
      
      write.csv(res_df, file = file.path(outdir, paste0("contrast_", cname, ".csv")), row.names = F)
      msg_ok(paste0("Saved: ", file.path(outdir, paste0("contrast_", cname, ".csv"))))
      
      assign(paste0("results","_contrast_", cname), res_df)
      
    }
    
    pca <- prcomp(t(assay(vsd)),scale. = T)
    pca_df <- as.data.frame(pca$x)
    pca_df$condition <- condition
    
    base_colors <- c("#804A45", "#455F80","#F5A553", "#BAC6D4", "#2E8B57", "#323840", "#FFF1C2")
    palette_auto <- base_colors[seq_len(length(unique(pca_df$condition)))]
    
    
    p <- ggplot(pca_df, aes(x=PC1, y=PC2, color=condition)) +
      scale_color_manual(values = palette_auto) +
      geom_point(size=3) +
      stat_ellipse(aes(group = condition), type = "norm", level = 0.95, linetype = "dashed") +
      labs(title="PCA",
           x=paste0("PC1 (", round(100*summary(pca)$importance[2,1], 1), "% var)"),
           y=paste0("PC2 (", round(100*summary(pca)$importance[2,2], 1), "% var)")) +
      theme_minimal() +
      geom_point(size=3) +
      stat_ellipse(aes(group = condition), type = "norm", level = 0.95, linetype = "dashed") +
      labs(title="PCA",
           x=paste0("PC1 (", round(100*summary(pca)$importance[2,1], 1), "% var)"),
           y=paste0("PC2 (", round(100*summary(pca)$importance[2,2], 1), "% var)")) +
      theme_minimal()
    
    ggsave(paste0(DEA_results_DIR,"/pca.png"), plot = p, width = 8, height = 6, dpi = 300)
    ggsave(paste0(DEA_results_DIR,"/pca.pdf"), plot = p, width = 8, height = 6, dpi = 300)
    
  }
}

## Generate Matrix for Prediction Model

if (method == "DESeq2") {
  expression_matrix <- assay(vsd)
  rownames(expression_matrix) <- gsub("#.*$", "", rownames(expression_matrix))
  
} else {
  mat.dge.log <- log2(mat.dge)
  rownames(mat.dge.log) <- gsub("#.*$", "", rownames(mat.dge.log))
  expression_matrix <- mat.dge.log}

write.csv(expression_matrix, paste0(DEA_results_DIR,"/expression_matrix.csv"))


###############################################################################
####                                                                       ####
####                VolcanoPlot & HeatMap & Barplot UP/DOWN                ####
####                                                                       ####
###############################################################################

res_names_filtered <- ls(pattern = "^res_filtered_contrast_")
empty_filtered <- res_names_filtered[sapply(res_names_filtered, function(x) {
  obj <- get(x)
  is.data.frame(obj) && nrow(obj) == 0
})]
empty_raw <- sub("^res_filtered_", "results_", empty_filtered)
res_names <- ls(pattern = "^results_contrast_")
res_names <- setdiff(res_names, empty_raw)

#res_names_filtered <- ls(pattern = "^res_filtered_contrast_")
res_names_filtered <- setdiff(res_names_filtered, empty_filtered)


if (length(res_names_filtered)==0){
  stop("None of the contrasts showed any differentially expressed features.")
}

# res_names es un vector de caracteres con nombres de data.frames
keep <- vapply(
  res_names_filtered,
  FUN = function(nm) nrow(get(nm)) >= 2,
  FUN.VALUE = logical(1)
)

res_names_filtered <- res_names_filtered[keep]

# raíz común: desde "contrast_" en adelante
root_res      <- sub("^results_",       "", res_names)
root_filtered <- sub("^res_filtered_",  "", res_names_filtered)

# quedarte con los res_names cuyo "root" también está en root_filtered
res_names_kept <- res_names[root_res %in% root_filtered]

# si quieres sobrescribir:
res_names <- res_names_kept
for (nm in res_names) {
  
message("Processing: ", nm)
  
results <- get(nm)

conds <- strsplit(sub("^results_contrast_", "", nm), "_vs_")[[1]]

## ---------------- VOLCANO PLOTS ----------------

if (method == "DESeq2"){
  log2FC <- results$log2FoldChange
  FDR <- results$padj
  dep.labels <- ifelse(log2FC > log2FC_thr & FDR < FDR_thr, paste0("Upregulated in ", unique(conds)[1]),
                       ifelse(log2FC < -log2FC_thr & FDR < FDR_thr, paste0("Downregulated in ", unique(conds)[1]), 
                              "Not significant"))
  
  volcano.df <- data.frame(
    RepeatSequence = rownames(results),
    log2FC = results$log2FoldChange,
    negLog10P = -log10(results$padj),
    DEG.Status = factor(dep.labels,
                        levels = c(paste0("Upregulated in ", unique(conds)[1]), 
                                   paste0("Downregulated in ", unique(conds)[1]),
                                   "Not significant")))
  
} else{
  log2FC <- results$logFC
  FDR <- results$FDR
  dep.labels <- ifelse(log2FC > log2FC_thr & FDR < FDR_thr, paste0("Upregulated in ", unique(conds)[1]),
                       ifelse(log2FC < -log2FC_thr & FDR < FDR_thr, paste0("Downregulated in ", unique(conds)[1]), 
                              "Not significant"))
  
  volcano.df <- data.frame(
    RepeatSequence = rownames(results),
    log2FC = results$logFC,
    negLog10P = -log10(results$FDR),
    DEG.Status = factor(dep.labels,
                        levels = c(paste0("Upregulated in ", unique(conds)[1]), 
                                   paste0("Downregulated in ", unique(conds)[1]),
                                   "Not significant")))
} 

## Check if there are DEGs
all_not_significant <- all(
  volcano.df$DEG.Status == "Not significant" | is.na(volcano.df$DEG.Status)
)

if (all_not_significant) {
  msg_error("No differentially expressed genes were detected. The pipeline will stop here.")
  next
}

grp_up   <- paste0("Upregulated in ", unique(conds)[1])
grp_down <- paste0("Downregulated in ", unique(conds)[1])
pal <- setNames(
  c("#804A45", "#455F80", "grey70"),          
  c(grp_up, grp_down, "Not significant")      
)

top_up   <- volcano.df |>  filter(DEG.Status == grp_up)    |>  slice_max(negLog10P, n = 10)
top_down <- volcano.df |> filter(DEG.Status == grp_down)  |>  slice_max(negLog10P, n = 10)
top_labs <- bind_rows(top_up, top_down)

write.table(top_labs, paste0(CWD,"Results/toplabels.txt"))
write.table(
  top_up$RepeatSequence,
  file = paste0(CWD,"Results/primersearchtoplabels.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


ggplot(volcano.df, aes(x = log2FC, y = negLog10P)) +
  geom_point(aes(color = DEG.Status), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = pal, breaks = names(pal))+
  geom_vline(xintercept = c(-log2FC_thr, log2FC_thr), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(FDR_thr), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_labs,
                  aes(label = RepeatSequence),
                  size = 2.5, max.overlaps = Inf, box.padding = 0.5, segment.size = 0.2) +
  labs(title = paste0("Volcano Plot (", unique(conds)[1], " vs ", unique(conds)[2], ")"),
       x = expression(log[2]*"(Fold-change)"),
       y = expression(-log[10]*"(padj)"),
       color = "DEG Status") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

ggsave(paste0(DEA_results_DIR,"/VolcanoPlot_", nm, ".png"),
       width = 6000, height = 4500, dpi = 600, units = "px")

ggsave(paste0(DEA_results_DIR,"/VolcanoPlot_", nm, ".pdf"),
       width = 6000, height = 4500, dpi = 600, units = "px")

## ---------------- HEATMAP ----------------

repeat_differentials <- subset(volcano.df, DEG.Status != "Not significant")  
repeat_differentials <- repeat_differentials[!is.na(repeat_differentials$DEG.Status),]
repeat_differentials <- repeat_differentials$RepeatSequence
repeat_differentials <- gsub("#.*$","", repeat_differentials)

expression_differentials <- expression_matrix[rownames(expression_matrix) %in% 
                                                repeat_differentials,]

annotation_col <- data.frame(
  Condition = factor(condition)
)
rownames(annotation_col) <- colnames(expression_differentials)

expression_differentials <- expression_differentials[,order(rownames(annotation_col))]

cond_levels <- levels(annotation_col$Condition)
base_colors <- c("#804A45", "#BAC6D4","#F5A553", "#2E8B57","#455F80", "#323840", "#FFF1C2")
ann_colors <- list(
  Condition = setNames(base_colors[seq_along(cond_levels)], cond_levels)
)

heat_colors <- colorRampPalette(c("#455F80", "white", "#804A45"))(100)

if (length(repeat_differentials) > 1){

heatmap <- pheatmap(expression_differentials,
                    annotation_col = annotation_col,
                    clustering_distance_rows = "euclidean",
                    cluster_cols = FALSE,
                    clustering_method = "complete",
                    scale = "row",
                    fontsize_row = if(length(expression_differentials) <= 30){8}else{4},
                    fontsize_col = if(length(samples) <= 16){8}else{6},
                    main = paste0("DEGs between conditions (", unique(conds)[1], " vs ", unique(conds)[2], ")"),
                    annotation_colors = ann_colors,
                    color = heat_colors,
                    border_color = NA,
                    annotation_names_col = if(length(samples) <= 30){TRUE}else{FALSE})
           }

png(paste0(DEA_results_DIR,"/heatmap_DEGs_", nm,".png"), width = 2600, height = 2000, res = 300)
grid.newpage()
grid.draw(heatmap$gtable)
dev.off()

pdf(paste0(DEA_results_DIR,"/heatmap_DEGs_", nm,".pdf"), width = 8.5, height = 6.5)
grid.newpage()
grid.draw(heatmap$gtable)
dev.off()

## ---------------- BARPLOT UP/DOWN ----------------

counts_deg <- volcano.df |> 
  filter(DEG.Status %in% c(grp_up, grp_down))

df_bar <- counts_deg %>%
  arrange(desc(log2FC)) %>%
  mutate(RepeatSequence = factor(RepeatSequence, levels = RepeatSequence))

ggplot(df_bar, aes(x = RepeatSequence, y = log2FC, fill = DEG.Status)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = pal) +
  labs(
    title = paste0("Differential Expression of Repetitive Elements (", unique(conds)[1], " vs ", unique(conds)[2],")"),
    x = "Repeat Sequence",
    y = "log2 Fold Change",
    fill = "Status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = if(length(counts_deg$RepeatSequence) <= 50){8}else{4}),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")

ggsave(paste0(DEA_results_DIR,"/BarPlotUpDown_", nm,".png"),
       width = 6000, height = 4500, dpi = 600, units = "px")

ggsave(paste0(DEA_results_DIR,"/BarPlotUpDown_", nm,".pdf"),
       width = 6000, height = 4500, dpi = 600, units = "px")

}



###############################################################################
####                                                                       ####
####                         Classification rep-types                      ####
####                                                                       ####
###############################################################################

## Read & Preprocess GTF
msg_info(paste0("Loading: ", repeatmasker_annotation_gtf))

gtf <- rtracklayer::readGFF(repeatmasker_annotation_gtf)

gtf_differentials <- gtf[gtf$gene_id %in% repeat_differentials,]
gtf_differentials <- as.data.frame(gtf_differentials)

repeat_class_info <- gtf_differentials |> 
  dplyr::select(gene_id, repeat_class) |> 
  distinct()

msg_info("Classifying repeats...")

for (nm in res_names) {

  
  results <- get(nm)
  
  conds <- strsplit(sub("^results_contrast_", "", nm), "_vs_")[[1]]
  
  if (method == "DESeq2"){
    log2FC <- results$log2FoldChange
    FDR <- results$padj
    dep.labels <- ifelse(log2FC > log2FC_thr & FDR < FDR_thr, paste0("Upregulated in ", unique(conds)[1]),
                         ifelse(log2FC < -log2FC_thr & FDR < FDR_thr, paste0("Downregulated in ", unique(conds)[1]), 
                                "Not significant"))
    
    volcano.df <- data.frame(
      RepeatSequence = rownames(results),
      log2FC = results$log2FoldChange,
      negLog10P = -log10(results$padj),
      DEG.Status = factor(dep.labels,
                          levels = c(paste0("Upregulated in ", unique(conds)[1]), 
                                     paste0("Downregulated in ", unique(conds)[1]),
                                     "Not significant")))
    
  } else{
    log2FC <- results$logFC
    FDR <- results$FDR
    dep.labels <- ifelse(log2FC > log2FC_thr & FDR < FDR_thr, paste0("Upregulated in ", unique(conds)[1]),
                         ifelse(log2FC < -log2FC_thr & FDR < FDR_thr, paste0("Downregulated in ", unique(conds)[1]), 
                                "Not significant"))
    
    volcano.df <- data.frame(
      RepeatSequence = rownames(results),
      log2FC = results$logFC,
      negLog10P = -log10(results$FDR),
      DEG.Status = factor(dep.labels,
                          levels = c(paste0("Upregulated in ", unique(conds)[1]), 
                                     paste0("Downregulated in ", unique(conds)[1]),
                                     "Not significant")))
  } 
  
  repeat_differentials <- subset(volcano.df, DEG.Status != "Not significant")
  repeat_differentials <- repeat_differentials$RepeatSequence
  repeat_differentials <- gsub("#.*$","", repeat_differentials)
  
  gtf_differentials <- gtf[gtf$gene_id %in% repeat_differentials,]
  gtf_differentials <- as.data.frame(gtf_differentials)
  
  repeat_class_info <- gtf_differentials |> 
    dplyr::select(gene_id, repeat_class) |> 
    distinct()  
  
DEGs_type <- volcano.df |>
  filter(DEG.Status != "Not significant") |> 
  dplyr::rename(gene_id = RepeatSequence) |> 
  mutate(gene_id = gsub("#.*$", "", gene_id)) |>    
  left_join(repeat_class_info, by = "gene_id")

type_df <- DEGs_type |>
  filter(!is.na(repeat_class)) |>
  dplyr::count(repeat_class) |>
  mutate(
    percentage = (n / sum(n)) * 100,
    percentage_label = paste0(round(percentage, 1), "%")
  ) |> 
  arrange(desc(n))


## Plot DEG
ggplot(type_df, aes(x = reorder(repeat_class, n), y = percentage)) +
  geom_bar(stat = "identity", fill = "#BAC6D4", alpha = 0.7) +
  geom_text(aes(label = percentage_label), 
            hjust = -0.1, size = 3.5) +
  coord_flip() +
  labs(
    title = paste0("Repeat class types distribution for ", unique(conds)[1], " vs ", unique(conds)[2]," - DEGs"),
    subtitle = paste("Total genes: ", sum(type_df$n)-length(which(duplicated(DEGs_type$gene_id)))),
    x = "Repeat class type",
    y = "Percentage (%)",
    caption = "DEA data analysis"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10)
  )

ggsave(paste0(DEA_results_DIR,"/Classification_DEGs_", nm, ".png"),
       width = 6000, height = 4500, dpi = 600, units = "px")

ggsave(paste0(DEA_results_DIR,"/Classification_DEGs_", nm, ".pdf"),
       width = 6000, height = 4500, dpi = 600, units = "px")

}


## Plot All
repeat_class_info_all <- gtf |> 
  dplyr::select(gene_id, repeat_class) |> 
  distinct()

rep_type <- volcano.df |> 
  dplyr::rename(gene_id = RepeatSequence) |> 
  mutate(gene_id = gsub("#.*$", "", gene_id)) |>    
  left_join(repeat_class_info_all, by = "gene_id")

rep_type <- rep_type |>
  dplyr::mutate(repeat_class = ifelse(is.na(repeat_class), "(Unknown)", repeat_class)) |>
  dplyr::count(repeat_class, name = "n") |>
  arrange(n) |>
  mutate(
    percentage = 100 * n / sum(n),
    percentage_label = sprintf("%.1f%%", percentage),
    repeat_class = factor(repeat_class, levels = repeat_class))

ggplot(rep_type, aes(x = repeat_class, y = percentage)) +
  geom_col(fill = "#BAC6D4", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = percentage_label),
            hjust = -0.1, size = 3.5) +
  coord_flip(clip = "off") +
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) +
  labs(
    title = paste0("Repeat class types distribution between all samples"),
    subtitle = paste("Total elements:", sum(rep_type$n)),
    x = "Repeat class type",
    y = "Percentage (%)",
    caption = "DEA data analysis"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.margin = margin(5.5, 30, 5.5, 5.5) 
  ) +
  expand_limits(y = max(rep_type$percentage) * 1.12)  

ggsave(paste0(DEA_results_DIR,"/Classification_All.png"),
       width = 6000, height = 4500, dpi = 600, units = "px")

ggsave(paste0(DEA_results_DIR,"/Classification_All.pdf"),
       width = 6000, height = 4500, dpi = 600, units = "px")

###############################################################################
####                                                                       ####
####                         Repeat RNAs distribution                      ####
####                                                                       ####
###############################################################################
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)})


if (method == "DESeq2"){log_means <- colMeans(assay(vsd), na.rm = TRUE)
} else{log_means <- log10(colMeans(mat.dge , na.rm = TRUE) + 1)} 

df_means <- data.frame(
  condition = sample_list$condition,
  mean_log = log_means
)

df_means$condition <- factor(df_means$condition)
my_comparisons <- combn(levels(df_means$condition), 2, simplify = FALSE)

base_colors <- c("#804A45", "#BAC6D4","#F5A553", "#2E8B57","#455F80", "#323840", "#FFF1C2")

ymax      <- max(df_means$mean_log)
y_bracket <- ymax * 1.08
y_label   <- ymax * 1.085

ggplot(df_means, aes(x = condition, y = mean_log, fill = condition)) +
  geom_violin(trim = FALSE, alpha = 0.6, width = 0.85, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8, color = "black") +
  geom_jitter(width = 0.12, size = 2, alpha = 0.9, color = "black") +
  scale_fill_manual(values = base_colors) +
  labs(
    title = "Mean repetitive sequence counts per replicate",
    x = "Condition",
    y = "Mean repetitive reads (log10)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  stat_compare_means(
    comparisons = my_comparisons,
    method      = "wilcox.test",
    label       = "p.signif")
  

ggsave(paste0(DEA_results_DIR,"/repetitive_counts_violin_box.png"), width = 8, height = 6, dpi = 300)
ggsave(paste0(DEA_results_DIR,"/repetitive_counts_violin_box.pdf"), width = 8, height = 6, dpi = 300)


## PDF Report

pdf_files <- list.files(DEA_results_DIR, pattern = "\\.pdf$", full.names = TRUE)

output_pdf <- file.path(DEA_results_DIR, "report.pdf")

if (length(pdf_files) > 0) {
  
  pdf_combine(input = pdf_files, output = output_pdf)
  msg_ok(paste0("PDF report has been generated successfully in ", DEA_results_DIR))
  
} else {
  msg_warn("No PDF files to create the report.")
}









