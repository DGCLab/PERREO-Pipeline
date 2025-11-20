args <- commandArgs(trailingOnly=TRUE)

GLOBAL_DIR <- args[[1]]
SAMPLES_DIR <- args[[2]]
REP_GTF_PATH <- args[[3]]
threads <- as.integer(args[[4]])
DEA_results <- args[[5]]
sample_list <- args[[6]]

print(GLOBAL_DIR)
print(REP_GTF_PATH)
print(threads)
print(DEA_results)


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
  message("✅ All the tables have the same rownames order.")
} else {
  message("⚠️ Tables do NOT have the same rownames order..")
  
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
  stop("❌ The condition is not met. Stopping the execution.")
}


