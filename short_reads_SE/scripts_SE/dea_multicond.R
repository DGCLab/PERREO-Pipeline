###############################################################################
####                                                                       ####
####                      DIFFERENTIAL EXPRESSION ANALYSIS                 ####
####                                                                       ####
####             Authors:                                                  ####
####                                                                       ####
####                                                                       ####
###############################################################################

# args <- commandArgs(trailingOnly=TRUE)
# 
# batch <- args[[1]]
# sample_list <- args[[2]]
# method <- args[[3]]
# CWD <- args[[4]]
# repeatmasker_annotation_gtf <- args[[5]]
# k_num <- args[[6]]
# 
# SAMPLES_DIR <- paste0(CWD,"/SAMPLES/")
# DEA_results_DIR <- paste0(SAMPLES_DIR,"/DEA_results")

method <- "DESeq2"
batch <- TRUE
k_num <- 2

#if (batch == "yes"){batch = TRUE}else{batch=FALSE}

## Load libraries
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

# Loading the metadata:

# sample_list <- read.table(paste0(CWD,"/",sample_list), header = T, sep="\t")
sample_list <- read.table("samplesheet_repeat_rna_3.txt", header = T, sep="\t")

## Prepare the data for DESeq2:

cts_raw <- read.table(paste0(SAMPLES_DIR,"count_data.txt"))
rownames <- rownames(cts_raw)

cts <- as.matrix(cts_raw)
rownames(cts) <- rownames
samples <- colnames(cts)

if (!identical(samples, sample_list$sample)) {
  stop("Order of samples ids in `cts` and `samplesheet` do not match!")
}

condition <- sample_list$condition
coldata <- data.frame(row.names = samples,
                      condition = factor(condition))

if (ncol(cts) != nrow(coldata)) {
  stop("Number of samples in `cts` and `coldata` do not match!")
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

outdir <- paste0(SAMPLES_DIR,"Contrasts")
dir.create(outdir, showWarnings = FALSE)

## Correct Batch Effect
if(batch == TRUE){
  ruvg <- TRUE
  
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
  
  set <- newSeqExpressionSet(as.matrix(mat),
                             phenoData = data.frame(coldata, row.names = colnames(cts)))
  set <- betweenLaneNormalization(set, which = "upper")  
  
  library(RUVSeq)
  res <- RUVg(x=set, cIdx=neg_controls, k=k_num ,isLog = F)
  mat.ruv <- res@assayData$normalizedCounts
  
  Wdf <- as.data.frame(pData(res)[, grep("^W_", colnames(pData(res))), drop = FALSE]) # Undesired factors (X columnas: W_1, W_2)
  
  # Add to metadata
  sample_list <- cbind(sample_list, Wdf)
  
  if (method == "edgeR") {
    library(edgeR)
    message("DEA is being performed by EdgeR")
    
    dge <- DGEList(counts = mat, samples = samples)
    dge <- calcNormFactors(dge, method="TMM")
    
    form <- as.formula(paste("~ ", paste(colnames(Wdf), collapse = " + "), "+ condition"))
    design <- model.matrix(form, data = sample_list)
    
    
    dge <- estimateDisp(dge, design)
    fit <- glmQLFit(dge, design)
    
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
        stop("Error: both conditions are baselines.")
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
      
      res_filtered <- res_df[abs(res_df$logFC) > 1 & res_df$FDR < 0.05, , drop = FALSE]
      assign(paste0("res_filtered","_contrast_", cname), res_filtered)
      
      #write.csv(res_df, file = file.path(outdir, paste0("contrast_", cname, ".csv")), row.names = F)
      # message("Saved: ", file.path(outdir, paste0("contrast_", cname, ".csv")))
      write.csv(res_df, file = paste0("contrast_", cname, ".csv"), row.names = F)
      message("Saved: ", paste0("contrast_", cname, ".csv"))
    }
  } else if (method == "DESeq2") {
    library(DESeq2)
    
    message("DEA is being performed by DESeq2")
    
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
      
      res_filtered <- subset(res_df, abs(log2FoldChange) > 1 & padj < 0.05) ## MIRAR
      assign(paste0("res_filtered","_contrast_", cname), res_filtered)
      
      # write.csv(res_df, file = file.path(outdir, paste0("contrast_", cname, ".csv")), row.names = F)
      # message("Saved: ", file.path(outdir, paste0("contrast_", cname, ".csv")))
      write.csv(res_df, file = paste0("contrast_", cname, ".csv"), row.names = F)
      message("Saved: ", paste0("contrast_", cname, ".csv"))
    }
    
    vsd <- varianceStabilizingTransformation(dds)
  }
  } else {
   if (method == "edgeR") {
    library(edgeR)
    message("DEA is being performed by EdgeR")
    
    dge <- DGEList(counts = mat, samples = samples)
    dge <- calcNormFactors(dge, method="TMM")
    
    design <- model.matrix(~ batch + condition, data = sample_list)
    dge <- estimateDisp(dge, design)
    fit <- glmQLFit(dge, design)
    
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
        stop("Error: both conditions are baselines.")
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
      
      res_filtered <- res_df[abs(res_df$logFC) > 1 & res_df$FDR < 0.05, , drop = FALSE]
      assign(paste0("res_filtered","_contrast_", cname), res_filtered)
      
      write.csv(res_df, file = file.path(outdir, paste0("contrast_", cname, ".csv")), row.names = F)
      message("Saved: ", file.path(outdir, paste0("contrast_", cname, ".csv")))
      
    }
  } else if (method == "DESeq2") 
    {
    library(DESeq2)
    
    message("DEA is being performed by DESeq2")
    
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
      
      res_filtered <- subset(res_df, abs(log2FoldChange) > 1 & padj < 0.05) ## MIRAR
      assign(paste0("res_filtered","_contrast_", cname), res_filtered)
      
      write.csv(res_df, file = file.path(outdir, paste0("contrast_", cname, ".csv")), row.names = F)
      message("Saved: ", file.path(outdir, paste0("contrast_", cname, ".csv")))
    }
    vsd <- varianceStabilizingTransformation(dds)
    }
  }
  
  #Once the batch effect is corrected we use the matrix for visualization
  #PCA
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
  
  
  pca_corrected <- prcomp(t(mat.ruv),scale. = T)
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

  } else {
  if (method == "edgeR") {
    library(edgeR)
    message("DEA is being performed by EdgeR")
    
    dge <- DGEList(counts = mat, samples = samples)
    dge <- calcNormFactors(dge, method="TMM")
    
    design <- model.matrix(~ condition, data = sample_list)
    dge <- estimateDisp(dge, design)
    fit <- glmQLFit(dge, design)
    
    mat.tmm <- cpm(mat)
    
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
        stop("Error: both conditions are baselines.")
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
      
      res_filtered <- res_df[abs(res_df$logFC) > 1 & res_df$FDR < 0.05, , drop = FALSE]
      assign(paste0("res_filtered","_contrast_", cname), res_filtered)
      
      write.csv(res_df, file = file.path(outdir, paste0("contrast_", cname, ".csv")), row.names = F)
      message("Saved: ", file.path(outdir, paste0("contrast_", cname, ".csv")))
      
      
      pca <- prcomp(t(assay(mat.tmm)),scale. = T)
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
    }
  } else if (method == "DESeq2") 
  {
    library(DESeq2)
    
    message("DEA is being performed by DESeq2")
    
    dds <- DESeqDataSetFromMatrix(
      countData = mat,
      colData = sample_list,
      design = ~ condition)
    
    dds <- DESeq(dds)
    
    for (p in pairs) {
      num <- p[2]   
      den <- p[1]   
      cname <- paste0(num, "_vs_", den)
      message("Processing contrast: ", cname)
      
      
      res <- results(dds, contrast = c("condition", num, den))
      res_df <- as.data.frame(res)
      res_df$gene <- rownames(res_df)
      
      res_filtered <- subset(res_df, abs(log2FoldChange) > 1 & padj < 0.05) ## MIRAR
      assign(paste0("res_filtered","_contrast_", cname), res_filtered)
      
      write.csv(res_df, file = file.path(outdir, paste0("contrast_", cname, ".csv")), row.names = F)
      message("Saved: ", file.path(outdir, paste0("contrast_", cname, ".csv")))
    }
    vsd <- varianceStabilizingTransformation(dds)
    
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
  }
}
