args <- commandArgs(trailingOnly=TRUE)

repeat_annotations <- args[[1]]
reference_genome <- args[[2]]
bam_dir <- args[[3]]
cwd <- args[[4]]

library(bambu)
library(ggbio)
library(ggplot2)
library(Rsamtools)
library(Biostrings)

genome <- readDNAStringSet(reference_genome)

names(genome) <- sub(" .*","",names(genome))

bambuAnnotations <- prepareAnnotations(repeat_annotations)

bam_paths <- list.files(path=bam_dir,pattern = "\\.bam$", full.names = TRUE, recursive=TRUE)

bamFiles <- BamFileList(bam_paths)

se.multisample <- bambu(reads = bamFiles, annotations = bambuAnnotations, 
                        genome = genome, stranded=TRUE)

saveRDS(se.multisample, file = paste0(cwd,"/se_multisample_bambu.rds"))

