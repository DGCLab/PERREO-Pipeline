args <- commandArgs(trailingOnly=TRUE)

GLOBAL_DIR <- args[[1]]
SAMPLES_DIR <- args[[2]]
REP_GTF_PATH <- args[[3]]
threads <- as.integer(args[[4]])
DEA_results <- args[[5]]
sample_list <- args[[6]]

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

# Required packages
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(DESeq2)

# Locate all the quant.txt inside Quantification folder 
files <- list.files(SAMPLES_DIR, pattern = "_?quant\\.txt$", recursive = TRUE, full.names = TRUE)
print(files)

global_matrix <- matrix()  

#Load the metadata file
file <- paste0(GLOBAL_DIR, "/", sample_list)
print(file)

meta <- read.table(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#Check that rownames are ordered in the same way in all the txt files
txt_files <- list()

for (i in 1:length(files)){
  txt <- read.table(files[i],sep="\t",header = T)
  txt_files[[i]] <- txt
  assign(colnames(txt),txt)
}

# Establish the order of the first matrix
orden_ref <- rownames(txt_files[[1]])

# Reorder the others
txt_files <- lapply(txt_files, function(tab) {
  tab[orden_ref, , drop = FALSE]
})

#Make sure if all of them follow the same order
todos_iguales <- all(sapply(txt_files, function(tab) {
  identical(rownames(tab), orden_ref)
}))

orden_ref <- rownames(txt_files[[1]])

todos_iguales <- all(sapply(txt_files, function(tab) {
  identical(rownames(tab), orden_ref)
}))

if (todos_iguales) {
  msg_ok("[FEATURECOUNTS] All the tables have the same rownames order.")
} else {
  msg_warn("[FEATURECOUNTS] Tables do NOT have the same rownames order..")
}


if (todos_iguales) {
  for (i in 1:length(files)){
    print(files[i])
    txt <- read.table(files[i],sep="\t",header = T)
    colnames(txt) <- sub("_.*", "", colnames(txt))
    global_matrix <- cbind(global_matrix,txt)
  }
  #Delete the original column of global_matrix
  global_matrix <- global_matrix[,-1]
  
  #Generate count matrix
  countData <- as.data.frame(global_matrix, check.names = FALSE)
  
  #Save the matrix
  write.table(countData, "count_data.txt",row.names = T)

}else{
  msg_error("[FEATURECOUNTS] The condition is not met. Stopping the execution.")
  stop("[FEATURECOUNTS] The condition is not met.", call. = FALSE)
}


