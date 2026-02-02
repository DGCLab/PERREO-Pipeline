args <- commandArgs(trailingOnly=TRUE)

MAP_DIR <- args[[1]]
sample_id <- args[[2]]
repeat_gtf <- args[[3]]
threads <- as.integer(args[[4]])
strandness <- args[[5]]
sample_dir <- args[[6]]
quant_dir <- args[[7]]

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


#####
# msg_info(MAP_DIR)
# msg_info(sample_id)
# msg_info(repeat_gtf)
# msg_info(quant_dir)

if (strandness=="unstranded"){
    strandness_fc = 0
}
if(strandness=="forward"){
  strandness_fc = 1
}
if(strandness=="reverse"){
  strandness_fc = 2
}

suppressPackageStartupMessages({
library(Rsubread)
library(openxlsx)
})

if (!file.exists(paste0(sample_dir,"_quant.txt"))){
     quant <- featureCounts(files = print(paste0(MAP_DIR,"/",sample_id,"_marked_duplicates_STAR.bam")), annot.ext = repeat_gtf,isGTFAnnotationFile = T, isPairedEnd = TRUE, GTF.attrType = "gene_id",countMultiMappingReads = TRUE,primaryOnly = FALSE, fraction=TRUE,strandSpecific = strandness_fc, nthreads=threads)
     msg_ok(paste0("[FEATURECOUNTS] ", sample_id,"quantification completed"))
     write.table(quant$counts,print(paste0(quant_dir,"/",sample_id,"_quant.txt")),sep="\t")
     write.table(quant$stat, print(paste0(quant_dir,"/",sample_id,"_quant_stats.txt")),sep="\t")
}
