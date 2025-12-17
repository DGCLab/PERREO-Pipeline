# Cargar el objeto featurecounts guardado
cts_raw <- read.table(paste0(SAMPLES_DIR,"/count_data.txt"))

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
  library(pdftools)  # Nuevo agregado para manejar PDFs.
})

# Resto del análisis DEA va aquí...

# === Consolidación de PDFs en un archivo de reporte ===
library(pdftools)

# Identificar todas las imágenes PDF en el directorio DEA_results_DIR
pdf_files <- list.files(DEA_results_DIR, pattern = "\\.pdf$", full.names = TRUE)

# Definir la ruta final del archivo reporte
output_pdf <- file.path(DEA_results_DIR, "report.pdf")

if (length(pdf_files) > 0) {
  # Combinar los archivos PDF en un único reporte
  pdf_combine(input = pdf_files, output = output_pdf)
  message(sprintf("Se ha generado el reporte en: %s", output_pdf))
} else {
  message("No se encontraron archivos PDF para consolidar en el reporte.")
}