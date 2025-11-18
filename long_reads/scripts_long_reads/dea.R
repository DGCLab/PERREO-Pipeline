###############################################################################
####                                                                       ####
####                      DIFFERENTIAL EXPRESSION ANALYSIS                 ####
####                                                                       ####
####             Authors:                                                  ####
####                                                                       ####
####                                                                       ####
###############################################################################

version <- "v1.0"
cat(sprintf("
  /^ ^\\
 / 0 0 \\              PERREO PIPELINE %s  
 V\\ Y /V      ----------------------------------
  / - \\       Processing Differential Expression
  |    \\                  Analysis
  || (__V)     
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

SAMPLES_DIR <- paste0(CWD,"/SAMPLES")
DEA_results_DIR <- paste0(SAMPLES_DIR,"/DEA_results")

# Cargar el objeto bambu guardado
se.multisample <- readRDS(paste0(CWD,"/se_multisample_bambu.rds"))

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
library(EDASeq)
})  

# Loading the metadata:

sample_list <- read.table(paste0(CWD,"/",sample_list), header = T, sep="\t")

## Prepare the data for DESeq2:

cts_raw <- assays(se.multisample)$counts
rownames <- rownames(cts_raw)

cts <- as.matrix(cts_raw)
rownames(cts) <- rownames
colnames(cts) <- gsub("\\.bam$", "", colnames(cts))
samples <- colnames(cts)

if (!identical(samples, sample_list$sample)) {
  sample_list <- sample_list[order(sample_list$sample),]
}

if (identical(as.character(samples), as.character(sample_list$sample))) {
  print("Order of samples ids in `cts` and `samplesheet` do match!")
} else {  
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

## Correct Batch Effect
if(batch == TRUE){
  
  ruvg <- TRUE
  message("Removing Batch Effect...")
  
  if ("batch" %in% colnames(sample_list)) {
    ruvg <- FALSE
  }
  
  if(ruvg == TRUE){
  #Extracting counts data to build the EdgeR object
  library(edgeR)
  message("Removing Batch Effect...")
  
  #Normalizing the object before applying RUVg algorithm
  
  design <- model.matrix(~condition)
  y <- DGEList(counts=mat, group=condition)
  y <- calcNormFactors(y, method="TMM")
  y <- estimateDisp(y, design)
  mat.tmm <- cpm(y)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef=colnames(design))
  
  #Selecting the 25% genes less statistically significant variables
  
  top_all <- topTags(lrt, n = nrow(mat))$table
  top_all <- top_all[order(top_all$PValue, decreasing = TRUE), ]
  
  #Selecting the 25% genes as stable genes
  
  n_empirical <- ceiling(0.25 * nrow(top_all))
  empirical <- rownames(top_all)[1:n_empirical]
  
  library(RUVSeq)
  res <- RUVg(x=log(mat.tmm+0.1), cIdx=empirical, k=k_num ,isLog = T)
  mat.tmm.ruv <- exp(res$normalizedCounts) ## Representation, WGCNA, Models...
  
  Wdf <- as.data.frame(pData(res)[, grep("^W_", colnames(pData(res))), drop = FALSE]) 
  
  # Add to metadata
  sample_list <- cbind(sample_list, Wdf)
  
  
  if (method == "edgeR") {
    library(edgeR)
    
    message("DEA is being performed by edgeR")
    
    # Crear DGEList con los conteos originales (no TMM)
    dge <- DGEList(counts = mat, samples = samples)
    
    # Diseño con factores RUV y condición
    form <- as.formula(paste("~ ", paste(colnames(Wdf), collapse = " + "), "+ condition"))
    design <- model.matrix(form, data = sample_list)
    
    dge <- calcNormFactors(dge, method="TMM")
    dge <- estimateDisp(dge, design)
    fit <- glmQLFit(dge, design)
    qlf <- glmQLFTest(fit, coef = colnames(design)[3])  # ajusta al nombre real
    
    results <- topTags(qlf, n = Inf)
    res_filtered <- results$table[abs(results$table$logFC)> log2FC_thr & results$table$FDR < FDR_thr,]
    
  } else if (method == "DESeq2") {
    library(DESeq2)
    
    message("DEA is being performed by DESeq2")
    
    form <- as.formula(paste("~ ", paste(colnames(Wdf), collapse = " + "), "+ condition"))
    
    dds <- DESeqDataSetFromMatrix(
      countData = mat,
      colData = sample_list,
      design = form
    )
    
    dds <- DESeq(dds)
    results <- results(dds, contrast = c("condition", unique(condition)[2], unique(condition)[1]))
    res_filtered <- subset(results, abs(log2FoldChange) > log2FC_thr & padj < FDR_thr) |> as.data.frame()
    
    vsd <- varianceStabilizingTransformation(dds)
  }
  } else {
    if (method == "edgeR") {
      library(edgeR)
      
      message("DEA is being performed by edgeR")
      
      dge <- DGEList(counts = mat, samples = samples)
      
      design <- model.matrix(~ batch + condition, data = sample_list)
      
      dge <- calcNormFactors(dge, method="TMM")
      dge <- estimateDisp(dge, design)
      fit <- glmQLFit(dge, design)
      qlf <- glmQLFTest(fit, coef = colnames(design)[3])  
      
      results <- topTags(qlf, n = Inf)
      res_filtered <- results$table[abs(results$table$logFC)> log2FC_thr & results$table$FDR < FDR_thr,]
      
    } else if (method == "DESeq2") {
      library(DESeq2)
      
      message("DEA is being performed by DESeq2")
      
      dds <- DESeqDataSetFromMatrix(
        countData = mat,
        colData = sample_list,
        design = ~ batch + condition
      )
      
      dds <- DESeq(dds)
      results <- results(dds, contrast = c("condition", unique(condition)[2], unique(condition)[1]))
      res_filtered <- subset(results, abs(log2FoldChange) > log2FC_thr & padj < FDR_thr) |> as.data.frame()
      
      vsd <- varianceStabilizingTransformation(dds)
    }
  }   
  #Once the batch effect is corrected we use the matrix for visualization
  #PCA
  pca <- prcomp(t(mat.tmm),scale. = T)
  pca_df <- as.data.frame(pca$x)
  pca_df$condition <- condition
  p <- ggplot(pca_df, aes(x=PC1, y=PC2, color=condition)) +
    geom_point(size=3) +
    scale_color_manual(values = c("#804A45", "#455F80")) +
    stat_ellipse(aes(group = condition), type = "norm", level = 0.95, linetype = "dashed") +
    labs(title="PCA - Before Batch Effect Corrected",
         x=paste0("PC1 (", round(100*summary(pca)$importance[2,1], 1), "% var)"),
         y=paste0("PC2 (", round(100*summary(pca)$importance[2,2], 1), "% var)")) +
    theme_minimal()
  
  ggsave(paste0(DEA_results_DIR,"pca_nobatch.png"), plot = p, width = 8, height = 6, dpi = 300)
  ggsave(paste0(DEA_results_DIR,"pca_nobatch.pdf"), plot = p, width = 8, height = 6, dpi = 300)
  
  pca_corrected <- prcomp(t(mat.tmm.ruv),scale. = T)
  pca_df_corrected <- as.data.frame(pca_corrected$x)
  pca_df_corrected$condition <- condition
  p_corrected <- ggplot(pca_df_corrected, aes(x=PC1, y=PC2, color=condition)) +
    geom_point(size=3) +
    scale_color_manual(values = c("#804A45", "#455F80")) +
    stat_ellipse(aes(group = condition), type = "norm", level = 0.95, linetype = "dashed") +
    labs(title="PCA - Batch Effect Corrected",
         x=paste0("PC1 (", round(100*summary(pca)$importance[2,1], 1), "% var)"),
         y=paste0("PC2 (", round(100*summary(pca)$importance[2,2], 1), "% var)")) +
    theme_minimal()
  
  ggsave(paste0(DEA_results_DIR,"/pca_corrected.png"), plot = p_corrected, width = 8, height = 6, dpi = 300)
  ggsave(paste0(DEA_results_DIR,"/pca_corrected.pdf"), plot = p_corrected, width = 8, height = 6, dpi = 300)
  
} else {
  if (method == "edgeR") {
    library(edgeR)
    
    message("DEA is being performed by edgeR")
    
    # Crear DGEList con los conteos originales (no TMM)
    dge <- DGEList(counts = mat, samples = samples)
    
    # Diseño con factores RUV y condición
    design <- model.matrix(~ condition, data = sample_list)
    
    dge <- calcNormFactors(dge, method="TMM")
    dge <- estimateDisp(dge, design)
    mat.tmm <- cpm(dge)
    fit <- glmQLFit(dge, design)
    qlf <- glmQLFTest(fit, coef = colnames(design)[2])  # ajusta al nombre real
    
    results <- topTags(qlf, n = Inf)
    res_filtered <- results$table[abs(results$table$logFC)> log2FC_thr & results$table$FDR < FDR_thr,]
    
    ## PCA
    
    pca <- prcomp(t(mat.tmm),scale. = T)
    pca_df <- as.data.frame(pca$x)
    pca_df$condition <- condition
    p <- ggplot(pca_df, aes(x=PC1, y=PC2, color=condition)) +
      geom_point(size=3) +
      scale_color_manual(values = c("#804A45", "#455F80")) +
      stat_ellipse(aes(group = condition), type = "norm", level = 0.95, linetype = "dashed") +
      labs(title="PCA",
           x=paste0("PC1 (", round(100*summary(pca)$importance[2,1], 1), "% var)"),
           y=paste0("PC2 (", round(100*summary(pca)$importance[2,2], 1), "% var)")) +
      theme_minimal()
    
    ggsave(paste0(DEA_results_DIR,"/pca.png"),plot = p, width = 8, height = 6, dpi = 300)
    ggsave(paste0(DEA_results_DIR,"/pca.pdf"),plot = p, width = 8, height = 6, dpi = 300)
    
  } else if (method == "DESeq2") {
    library(DESeq2)
    
    message("DEA is being performed by DESeq2")
    
    dds <- DESeqDataSetFromMatrix(
      countData = mat,
      colData = sample_list,
      design = ~ condition
    )
    
    dds <- DESeq(dds)
    results <- results(dds, contrast = c("condition", unique(condition)[2], unique(condition)[1]))
    res_filtered <- subset(results, abs(log2FoldChange) > log2FC_thr & padj < FDR_thr) |> as.data.frame()
    
    vsd <- varianceStabilizingTransformation(dds)
    
    ## PCA
    pca <- prcomp(t(assay(vsd)),scale. = T)
    pca_df <- as.data.frame(pca$x)
    pca_df$condition <- condition
    p <- ggplot(pca_df, aes(x=PC1, y=PC2, color=condition)) +
      scale_color_manual(values = c("#804A45", "#455F80")) +
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


if (batch == TRUE){
  mat.tmm.ruv.log <- log2(mat.tmm.ruv)
  rownames(mat.tmm.ruv.log) <- gsub("#.*$","",rownames(mat.tmm.ruv.log))
  expression_matrix <- mat.tmm.ruv.log
} else{
  if (method == "DESeq2") {
  vsd <- assay(vsd)
  rownames(vsd) <- gsub("#.*$","",rownames(vsd))
  expression_matrix <- vsd} else {
    mat.tmm.log <- log2(mat.tmm)
    rownames(mat.tmm.log) <- gsub("#.*$", "", rownames(mat.tmm.log))
    expression_matrix <- mat.tmm.log}
} 

write.csv(expression_matrix, paste0(DEA_results_DIR,"/expression_matrix.csv"))
write.csv(res_filtered, paste0(DEA_results_DIR,"/DEG.csv"))
write.csv(results, paste0(DEA_results_DIR,"/DEA_results.csv"))


###############################################################################
####                                                                       ####
####                VolcanoPlot & HeatMap & Barplot UP/DOWN                ####
####                                                                       ####
###############################################################################


## ---------------- VOLCANO PLOTS ----------------
results <- as.data.frame(results)

if (method == "DESeq2"){
  log2FC <- results$log2FoldChange
  FDR <- results$padj
  dep.labels <- ifelse(log2FC > log2FC_thr & FDR < FDR_thr, paste0("Upregulated in ", unique(condition)[2]),
                       ifelse(log2FC < -log2FC_thr & FDR < FDR_thr, paste0("Downregulated in ", unique(condition)[2]), 
                              "Not significant"))
  
  volcano.df <- data.frame(
    RepeatSequence = rownames(results),
    log2FC = results$log2FoldChange,
    negLog10P = -log10(results$padj),
    DEG.Status = factor(dep.labels,
                        levels = c(paste0("Upregulated in ", unique(condition)[2]), 
                                   paste0("Downregulated in ", unique(condition)[2]),
                                   "Not significant")))
  
} else{
  log2FC <- results$logFC
  FDR <- results$FDR
  dep.labels <- ifelse(log2FC > log2FC_thr & FDR < FDR_thr, paste0("Upregulated in ", unique(condition)[2]),
                       ifelse(log2FC < -log2FC_thr & FDR < FDR_thr, paste0("Downregulated in ", unique(condition)[2]), 
                              "Not significant"))
  
  volcano.df <- data.frame(
    RepeatSequence = rownames(results),
    log2FC = results$logFC,
    negLog10P = -log10(results$FDR),
    DEG.Status = factor(dep.labels,
                        levels = c(paste0("Upregulated in ", unique(condition)[2]), 
                                   paste0("Downregulated in ", unique(condition)[2]),
                                   "Not significant")))
} 

## Check if there are DEGs
all_not_significant <- all(
  volcano.df$DEG.Status == "Not significant" | is.na(volcano.df$DEG.Status)
)

if (all_not_significant) {
  message("No differentially expressed genes were detected. The pipeline will stop here.")
  quit(save = "no", status = 0)
}

grp_up   <- paste0("Upregulated in ", unique(condition)[2])
grp_down <- paste0("Downregulated in ", unique(condition)[2])
pal <- setNames(
  c("#804A45", "#455F80", "grey70"),          
  c(grp_up, grp_down, "Not significant")      
)

top_up   <- volcano.df |>  filter(DEG.Status == grp_up)    |>  slice_max(negLog10P, n = 10)
top_down <- volcano.df |> filter(DEG.Status == grp_down)  |>  slice_max(negLog10P, n = 10)
top_labs <- bind_rows(top_up, top_down)
write.table(top_labs, "toplabels.txt")

ggplot(volcano.df, aes(x = log2FC, y = negLog10P)) +
  geom_point(aes(color = DEG.Status), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = pal, breaks = names(pal))+
  geom_vline(xintercept = c(-log2FC_thr, log2FC_thr), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(FDR_thr), linetype = "dashed", color = "black") +
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

ggsave(paste0(DEA_results_DIR,"/VolcanoPlot.pdf"),
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
rownames(annotation_col) <- colnames(expression_differentials)


cond_levels <- levels(annotation_col$Condition)
base_colors <- c("#804A45", "#BAC6D4","#F5A553", "#2E8B57","#455F80", "#323840", "#FFF1C2")
ann_colors <- list(
  Condition = setNames(base_colors[seq_along(cond_levels)], cond_levels)
)

heat_colors <- colorRampPalette(c("#455F80", "white", "#804A45"))(100)

heatmap <- pheatmap(expression_differentials,
                    annotation_col = annotation_col,
                    clustering_distance_rows = "euclidean",
                    cluster_cols = FALSE,
                    clustering_method = "complete",
                    scale = "row",
                    fontsize_row = if(length(expression_differentials) <= 30){8}else{4},
                    fontsize_col = if(length(samples) <= 16){8}else{6},
                    main = "DEGs between conditions",
                    annotation_colors = ann_colors,
                    color = heat_colors,
                    border_color = NA,
                    annotation_names_col = if(length(samples) <= 30){TRUE}else{FALSE})

png(paste0(DEA_results_DIR,"/heatmap_DEGs.png"), width = 2600, height = 2000, res = 300)
grid.newpage()
grid.draw(heatmap$gtable)
dev.off()

pdf(paste0(DEA_results_DIR,"/heatmap_DEGs.pdf"), width = 8.5, height = 6.5)
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
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = if(length(counts_deg$RepeatSequence) <= 50){8}else{4}),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")

ggsave(paste0(DEA_results_DIR,"/BarPlotUpDown.png"),
       width = 6000, height = 4500, dpi = 600, units = "px")

ggsave(paste0(DEA_results_DIR,"/BarPlotUpDown.pdf"),
       width = 6000, height = 4500, dpi = 600, units = "px")

###############################################################################
####                                                                       ####
####                         Classification rep-types                      ####
####                                                                       ####
###############################################################################

## Read & Preprocess GTF
message("Loading: ", repeatmasker_annotation_gtf)

gtf <- rtracklayer::readGFF(repeatmasker_annotation_gtf)
gtf_df_clean <- gtf |> 
  mutate(repeat_family = stringr::str_extract(transcript_id, "(?<=#)[^/]+(?=/|$)"),
         gene_id = stringr::str_extract(gene_id, "^[^#]+"), 
         repeat_family = stringr::str_remove(repeat_family, "_dup.*"))

gtf_differentials <- gtf_df_clean[gtf_df_clean$gene_id %in% repeat_differentials,]
gtf_differentials <- as.data.frame(gtf_differentials)

message("Classificating repeats...")

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

ggsave(paste0(DEA_results_DIR,"/Classification_DEGs.pdf"),
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

ggsave(paste0(DEA_results_DIR,"/Classification_All.pdf"),
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
ggsave(paste0(DEA_results_DIR,"/repetitive_counts_violin_box.pdf"), plot = p, width = 8, height = 6, dpi = 300)


