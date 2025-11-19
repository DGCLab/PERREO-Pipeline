args <- commandArgs(trailingOnly=TRUE)

DEA_DIR <- args[[1]]
cwd <- args[[2]]
sample_list <- args[[3]]
coexpression_dir <- args[[4]]

library(WGCNA)
library(tidyverse)     
library(magrittr) 
library(scales)
library(stringr)

samplesheet <- read.table(paste0(cwd,"/",sample_list),sep="\t",header = T)

expression_matrix <- read.csv2(paste0(DEA_DIR,"/expression_matrix.csv"),sep=",",header=T)
expression_matrix <- as.data.frame(expression_matrix)
expression_matrix <-expression_matrix[!duplicated(expression_matrix$X),]
rownames(expression_matrix) <- expression_matrix$X
expression_matrix <- expression_matrix[,-1]
expression_matrix <- as.matrix(expression_matrix)

expression = t(expression_matrix)

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(expression, powerVector = powers, verbose = 5)

#Plotting the results
#sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9


targetR2 <- 0.9
softPower <- sft$fitIndices$Power[which(sft$fitIndices$SFT.R.sq > targetR2)[1]]
if (is.na(softPower)) {
  softPower <- max(sft$fitIndices$Power)
  message("⚠️ Ninguna potencia alcanza R² ≥ ", targetR2, ". Usando la máxima: ", softPower)
} else {
  message("✅ SoftPower elegido automáticamente: ", softPower)
}

adjacency = adjacency(expression, power = softPower, type = "unsigned") #Calculating the adjacency matrix
#help(adjacency )

#TOM = TOMsimilarity(adjacency) #Calculating the topological overlap matrix

#dissTOM = 1-TOM ##Calculating the dissimilarity



temp_cor <- cor       
cor <- WGCNA::cor 

storage.mode(expression) <- "numeric"

netwk <- blockwiseModules(expression,               
                          power = softPower,                
                          networkType = "signed",
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          numericLabels = T,
                          verbose = 3)

# For plotting it is necessary to convert to colors
pdf(paste0(coexpression_dir,"/modules_dendrogram.pdf"), width = 8, height = 8)

mergedColors = labels2colors(netwk$colors)
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

dev.off()

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

write_delim(module_df,
            file = paste0(coexpression_dir,"/gene_modules.txt"),
            delim = "\t")


# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(expression,mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$sample = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-sample) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME <- mME %>%
  mutate(sample = factor(sample))

# We keep one row per sample
samplesheet_simple <- samplesheet %>%
  distinct(sample, condition)

cond_levels <- samplesheet_simple %>%
  pull(condition) %>%
  unique()

cond_palette <- hue_pal()(length(cond_levels))

cond_cols <- setNames(cond_palette, cond_levels)

cond_por_sample <- samplesheet_simple$condition[
  match(levels(mME$sample), samplesheet_simple$sample)
]

axis_cols <- cond_cols[cond_por_sample]


pdf(paste0(coexpression_dir,"/modules_samples_heatmap.pdf"), width = 8, height = 8)

mME %>%
  ggplot(aes(x = sample, y = name, fill = value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1, 1)
  ) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      color = axis_cols   
    )
  ) +
  labs(
    title = "Module-sample Relationships",
    y = "Modules",
    fill = "corr"
  )

dev.off()

ss <- samplesheet %>%
  distinct(sample, condition)

# Reordering ss in order to keep the same order as MEs0
ss <- ss[match(rownames(MEs0), ss$sample), ]

# Basic check
stopifnot(all(rownames(MEs0) == ss$sample))

ME_cols <- grep("^ME", colnames(MEs0), value = TRUE)
MEs_mat <- MEs0[, ME_cols, drop = FALSE]

cond_factor <- factor(ss$condition)
cond_levels <- levels(cond_factor)

datTraits <- sapply(cond_levels, function(cl) as.numeric(cond_factor == cl))
datTraits <- as.data.frame(datTraits)
colnames(datTraits) <- cond_levels
rownames(datTraits) <- ss$sample  

#Correlations and p-values
moduleTraitCor  <- cor(MEs_mat, datTraits, use = "p")
moduleTraitPval <- corPvalueStudent(moduleTraitCor, nrow(MEs_mat))

textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPval, 1), ")")

#Heatmap of correlation between conditions and modules
pdf(paste0(coexpression_dir,"/labeled_heatmap.pdf"), width = 8, height = 8)

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits),      
               yLabels = colnames(MEs_mat),        
               ySymbols = colnames(MEs_mat),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = TRUE,
               cex.text = 0.8,
               zlim = c(-1, 1),
               main = "Modules and experimental conditions correlation")

dev.off()


#Generation of networks of interest

# Selecting the top correlation values
moduleMaxCor <- apply(moduleTraitCor, 1, function(x) max(abs(x), na.rm = TRUE))

# Ordering from the highest to the lowest correlation value
moduleMaxCor <- sort(moduleMaxCor, decreasing = TRUE)

# Keeping the top 3
nTop <- min(3, length(moduleMaxCor))
topModulesME <- names(moduleMaxCor)[1:nTop]   # p.ej. "MEturquoise", "MEblue", ...

topModuleColors <- gsub("^ME", "", topModulesME)

topModulesME
topModuleColors

mergedColors = labels2colors(netwk$colors)
# names(mergedColors) deberían ser los IDs de los genes, mismo orden que en expression
geneNames <- names(mergedColors)

library(WGCNA)

# Asegúrate de que adjacency y mergedColors están en el mismo orden de genes
stopifnot(all(rownames(adjacency) == geneNames),
          all(colnames(adjacency) == geneNames))

for (col in topModuleColors) {
  message("Processing module: ", col)
  
  # Genes pertenecientes a este módulo
  inModule <- mergedColors == col
  modGenes <- geneNames[inModule]
  
  # Submatriz de adyacencia para el módulo
  modAdj <- adjacency[inModule, inModule]
  
  # Opcional: umbral de adyacencia para la exportación (ajusta a tu gusto)
  thr <- 0.1
  
  # Exportar a Cytoscape
  exportNetworkToCytoscape(
    modAdj,
    edgeFile = paste0(coexpression_dir,"/CytoscapeInput-edges-", col, ".txt"),
    nodeFile = paste0(coexpression_dir,"/CytoscapeInput-nodes-", col, ".txt"),
    weighted = TRUE,
    threshold = thr,
    nodeNames = modGenes,
    altNodeNames = modGenes,
    nodeAttr = mergedColors[inModule]
  )
}



