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


# Paquetes mínimos
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(DESeq2)

# 1) Localiza todos los quant.txt dentro de */Quantification/
files <- list.files(SAMPLES_DIR, pattern = "_?quant\\.txt$", recursive = TRUE, full.names = TRUE)
print(files)

global_matrix <- matrix()  # 3 filas, 0 columnas


#Load the metadata file
file <- paste0(GLOBAL_DIR, "/", sample_list)
print(file)
#file.exists(file)

meta <- read.table(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#Primero compruebo que los rownames son idénticos en orden en todos los txt
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
  message("✅ Todas las tablas tienen el mismo orden de rownames.")
} else {
  message("⚠️ Las tablas NO tienen el mismo orden de rownames.")
  
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
  
  # 5) Genera la matriz de conteo
  countData <- as.data.frame(global_matrix, check.names = FALSE)
  
  #(opcional) guarda la matriz y el objeto
  write.table(countData, paste0(GLOBAL_DIR,"/Results/count_data.txt"),row.names = T)

}else{
  stop("❌ La condición no se cumple. Deteniendo la ejecución.")
}


