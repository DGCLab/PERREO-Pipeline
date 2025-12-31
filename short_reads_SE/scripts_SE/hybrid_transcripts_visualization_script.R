args <- commandArgs(trailingOnly=TRUE)

CWD <- args[[1]]
sample_list <- args[[2]]

suppressPackageStartupMessages({library(ggpubr)})

hybrid_transcripts <- read.table(paste0(CWD,"/Results/transcriptome_assembly/hybrid_transcripts_summary.tsv"),sep="\t",header=T)
sample_list <- read.table(paste0(CWD,"/",sample_list),sep="\t",header=T)

hybrid_transcripts$sample <- gsub("_transcriptome","",hybrid_transcripts$sample)

hybrid_transcripts <- merge(hybrid_transcripts, sample_list,
                     by = "sample")


hybrid_transcripts <- hybrid_transcripts[order(hybrid_transcripts$condition),]
hybrid_transcripts$hybrid_frac <- gsub(",",".",hybrid_transcripts$hybrid_frac)
hybrid_transcripts$hybrid_frac <- as.numeric(hybrid_transcripts$hybrid_frac)


suppressPackageStartupMessages({library(ggplot2)})

base_colors <- c("#804A45", "#BAC6D4","#F5A553", "#2E8B57","#455F80", "#323840", "#FFF1C2")

my_comparisons <- combn(unique(sample_list$condition), 2, simplify = FALSE)


p <- ggplot(hybrid_transcripts, aes(x = condition, y = hybrid_frac, fill = condition)) +
  geom_violin(trim = FALSE, alpha = 0.6) +          
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +       
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) + 
  stat_summary(fun = mean, geom = "point",          
               color = "black", size = 3) +
  scale_fill_manual(values = base_colors) +
  theme_bw(base_size = 12) +
  labs(x = "Condition", y = "Percentage of hybrid transcripts") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  stat_compare_means(
    comparisons = my_comparisons,
    method      = "wilcox.test",
    label       = "p.signif")

ggsave(paste0(CWD,"/Results/transcriptome_assembly/hybrid_transcripts.pdf"), plot = p, width = 8, height = 6, dpi = 300)

