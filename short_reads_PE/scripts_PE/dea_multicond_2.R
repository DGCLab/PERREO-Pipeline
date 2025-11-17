###############################################################################
####                                                                       ####
####                      DIFFERENTIAL EXPRESSION ANALYSIS                 ####
####                                                                       ####
####             Authors:                                                  ####
####                                                                       ####
####                                                                       ####
###############################################################################

args <- commandArgs(trailingOnly=TRUE)

batch <- args[[1]]
sample_list <- args[[2]]
method <- args[[3]]
CWD <- args[[4]]
repeatmasker_annotation_gtf <- args[[5]]
k_num <- args[[6]]
FDR <- args[[7]]
log2FC <- args[[8]]
 
SAMPLES_DIR <- paste0(CWD,"/SAMPLES/")
DEA_results_DIR <- paste0(SAMPLES_DIR,"/DEA_results")

#method <- "DESeq2"
#batch <- TRUE
#k_num <- 2

if (batch == "yes"){batch = TRUE}else{batch=FALSE}

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


sample_list <- read.table(paste0(CWD,"/",sample_list), header = T, sep="\t")

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
      
      res_filtered <- res_df[abs(res_df$logFC) > "$log2FC" & res_df$FDR < "$FDR", , drop = FALSE]
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
      
      res_filtered <- subset(res_df, abs(log2FoldChange) > "$log2FC" & padj < "$FDR") ## MIRAR
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
      
      res_filtered <- res_df[abs(res_df$logFC) > "$log2FC" & res_df$FDR < "$FDR", , drop = FALSE]
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
      
      res_filtered <- subset(res_df, abs(log2FoldChange) > "$log2FC" & padj < "$FDR") ## MIRAR
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
      
      res_filtered <- res_df[abs(res_df$logFC) > "$log2FC" & res_df$FDR < "$FDR", , drop = FALSE]
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
      
      res_filtered <- subset(res_df, abs(log2FoldChange) > "$log2FC" & padj < "$FDR") ## MIRAR
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


## ---------------- VOLCANO PLOTS ----------------
results <- as.data.frame(results)

if (method == "DESeq2"){
  log2FC <- results$log2FoldChange
  FDR <- results$padj
  dep.labels <- ifelse(log2FC > 1 & FDR < 0.05, paste0("Upregulated in ", unique(condition)[1]),
                       ifelse(log2FC < -1 & FDR < 0.05, paste0("Downregulated in ", unique(condition)[1]), 
                              "Not significant"))
  
  volcano.df <- data.frame(
    RepeatSequence = rownames(results),
    log2FC = results$log2FoldChange,
    negLog10P = -log10(results$padj),
    DEG.Status = factor(dep.labels,
                        levels = c(paste0("Upregulated in ", unique(condition)[1]), 
                                   paste0("Downregulated in ", unique(condition)[1]),
                                   "Not significant")))
  
} else{
  log2FC <- results$logFC
  FDR <- results$FDR
  dep.labels <- ifelse(log2FC > 1 & FDR < 0.05, paste0("Upregulated in ", unique(condition)[1]),
                       ifelse(log2FC < -1 & FDR < 0.05, paste0("Downregulated in ", unique(condition)[1]), 
                              "Not significant"))
  
  volcano.df <- data.frame(
    RepeatSequence = rownames(results),
    log2FC = results$logFC,
    negLog10P = -log10(results$FDR),
    DEG.Status = factor(dep.labels,
                        levels = c(paste0("Upregulated in ", unique(condition)[1]), 
                                   paste0("Downregulated in ", unique(condition)[1]),
                                   "Not significant")))
} 


grp_up   <- paste0("Upregulated in ", unique(condition)[1])
grp_down <- paste0("Downregulated in ", unique(condition)[1])
pal <- setNames(
  c("#804A45", "#455F80", "grey70"),          
  c(grp_up, grp_down, "Not significant")      
)

top_up   <- volcano.df |>  filter(DEG.Status == grp_up)    |>  slice_max(negLog10P, n = 10)
top_down <- volcano.df |> filter(DEG.Status == grp_down)  |>  slice_max(negLog10P, n = 10)
top_labs <- bind_rows(top_up, top_down)

ggplot(volcano.df, aes(x = log2FC, y = negLog10P)) +
  geom_point(aes(color = DEG.Status), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = pal, breaks = names(pal))+
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_labs,
                  aes(label = RepeatSequence),
                  size = 2.5, max.overlaps = Inf, box.padding = 0.5, segment.size = 0.2) +
  labs(title = paste0("Volcano Plot (", unique(condition)[1], " vs ", unique(condition)[2], ")"),
       x = expression(log[2]*"(Fold-change)"),
       y = expression(-log[10]*"(padj)"),
       color = "DEG Status") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

ggsave(paste0(DEA_results_DIR,"/VolcanoPlot.png"),
       width = 6000, height = 4500, dpi = 600, units = "px")


## ---------------- HEATMAP ----------------

repeat_differentials <- subset(volcano.df, DEG.Status != "Not significant")
repeat_differentials <- repeat_differentials$RepeatSequence
repeat_differentials <- gsub("#.*$","", repeat_differentials)

expression_differentials <- expression_matrix[rownames(expression_matrix) %in% 
                                                repeat_differentials,]


annotation_col <- data.frame(
  Condition = factor(condition, levels = c(unique(condition)[2], unique(condition)[1]))
)
rownames(annotation_col) <- colnames(mat)

ann_colors <- list(
  Condition = c(ESO = "#804A45", HC = "#BAC6D4")
)

heat_colors <- colorRampPalette(c("#804A45", "white", "#455F80"))(100)

heatmap <- pheatmap(expression_differentials,
                    annotation_col = annotation_col,
                    clustering_distance_rows = "euclidean",
                    cluster_cols = FALSE,
                    clustering_method = "complete",
                    scale = "row",
                    fontsize_row = 8,
                    fontsize_col = 10,
                    main = "DEGs between conditions",
                    annotation_colors = ann_colors,
                    color = heat_colors)

png(paste0(DEA_results_DIR,"/heatmap_DEGs.png"), width = 2600, height = 2000, res = 300)
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
    title = "Differential Expression of Repetitive Elements",
    x = "Repeat Sequence",
    y = "log2 Fold Change",
    fill = "Status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")

ggsave(paste0(DEA_results_DIR,"/BarPlotUpDown.png"),
       width = 6000, height = 4500, dpi = 600, units = "px")


###############################################################################
####                                                                       ####
####                         Classification rep-types                      ####
####                                                                       ####
###############################################################################

## Read & Preprocess GTF
gtf <- rtracklayer::readGFF("repeatmasker_nanopore_annotation.gtf")
gtf_df_clean <- gtf |> 
  mutate(repeat_family = stringr::str_extract(transcript_id, "(?<=#)[^/]+(?=/|$)"),
         gene_id = stringr::str_extract(gene_id, "^[^#]+"), 
         repeat_family = stringr::str_remove(repeat_family, "_dup.*"))

gtf_differentials <- gtf_df_clean[gtf_df_clean$gene_id %in% repeat_differentials,]
gtf_differentials <- as.data.frame(gtf_differentials)

repeat_class_info <- gtf_differentials |> 
  dplyr::select(gene_id, repeat_family) |> 
  distinct()

DEGs_type <- volcano.df %>%
  filter(DEG.Status != "Not significant") |> 
  dplyr::rename(gene_id = RepeatSequence) |> 
  mutate(gene_id = gsub("#.*$", "", gene_id)) |>    
  left_join(repeat_class_info, by = "gene_id")

type_df <- DEGs_type %>%
  # Filtrar genes que tienen repeat_class asignado
  filter(!is.na(repeat_family)) %>%
  # Contar cada tipo
  dplyr::count(repeat_family) %>%
  # Calcular porcentajes
  mutate(
    percentage = (n / sum(n)) * 100,
    percentage_label = paste0(round(percentage, 1), "%")
  ) %>%
  # Ordenar por frecuencia (descendente)
  arrange(desc(n))


## Plot DEG
ggplot(type_df, aes(x = reorder(repeat_family, n), y = percentage)) +
  geom_bar(stat = "identity", fill = "#BAC6D4", alpha = 0.7) +
  geom_text(aes(label = percentage_label), 
            hjust = -0.1, size = 3.5) +
  coord_flip() +
  labs(
    title = "Repeat class types distribution - DEGs",
    subtitle = paste("Total genes: ", sum(type_df$n)),
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

ggsave(paste0(DEA_results_DIR,"/Classification_DEGs.png"),
       width = 6000, height = 4500, dpi = 600, units = "px")


## Plot All
repeat_class_info_all <- gtf_df_clean |> 
  dplyr::select(gene_id, repeat_family) |> 
  distinct()

rep_type <- volcano.df |> 
  dplyr::rename(gene_id = RepeatSequence) |> 
  mutate(gene_id = gsub("#.*$", "", gene_id)) |>    
  left_join(repeat_class_info_all, by = "gene_id")

rep_type <- rep_type %>%
  dplyr::mutate(repeat_family = ifelse(is.na(repeat_family), "(Unknown)", repeat_family)) %>%
  dplyr::count(repeat_family, name = "n") %>%
  arrange(n) %>%
  mutate(
    percentage = 100 * n / sum(n),
    percentage_label = sprintf("%.1f%%", percentage),
    repeat_family = factor(repeat_family, levels = repeat_family))

ggplot(rep_type, aes(x = repeat_family, y = percentage)) +
  geom_col(fill = "#BAC6D4", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = percentage_label),
            hjust = -0.1, size = 3.5) +
  coord_flip(clip = "off") +
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) +
  labs(
    title = "Repeat class types distribution",
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


###############################################################################
####                                                                       ####
####                         Repeat RNAs distribution                      ####
####                                                                       ####
###############################################################################

library(ggplot2)
library(ggpubr)

if (method == "DESeq2"){log_means <- colMeans(vsd, na.rm = TRUE)
} else{log_means <- log10(colMeans(mat.tmm , na.rm = TRUE) + 1)} 

df_means <- data.frame(
  condition = sample_list$condition,
  mean_log = log_means
)

ymax      <- max(df_means$mean_log)
y_bracket <- ymax * 1.08
y_label   <- ymax * 1.085

p <- ggplot(df_means, aes(x = condition, y = mean_log, fill = condition)) +
  geom_violin(trim = FALSE, alpha = 0.6, width = 0.85, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8, color = "black") +
  geom_jitter(width = 0.12, size = 2, alpha = 0.9, color = "black") +
  
  scale_fill_manual(values = c("#804A45", "#BAC6D4")) +
  labs(
    title = "Mean repetitive sequence counts per replicate",
    x = "Condition",
    y = "Mean repetitive reads (log10)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    label.x = 1.5,
    label.y = y_label
  ) +
  geom_segment(aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket), linewidth = 0.7) +
  geom_segment(aes(x = 1, xend = 1, y = y_bracket, yend = y_bracket - 0.02), linewidth = 0.7) +
  geom_segment(aes(x = 2, xend = 2, y = y_bracket, yend = y_bracket - 0.02), linewidth = 0.7)

ggsave(paste0(DEA_results_DIR,"/repetitive_counts_violin_box.png"), plot = p, width = 8, height = 6, dpi = 300)
