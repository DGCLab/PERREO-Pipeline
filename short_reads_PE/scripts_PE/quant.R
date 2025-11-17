args <- commandArgs(trailingOnly=TRUE)

MAP_DIR <- args[[1]]
sample_id <- args[[2]]
repeat_gtf <- args[[3]]
threads <- as.integer(args[[4]])
strandness <- args[[5]]
sample_dir <- args[[6]]
quant_dir <- args[[7]]

print(MAP_DIR)
print(sample_id)
print(repeat_gtf)
print(threads)
print(quant_dir)

if(strandness=="forward"){
  strandness_fc = 1
}
if(strandness=="reverse"){
  strandness_fc = 2
}

print(strandness)
print(strandness_fc)

library(Rsubread)
library(openxlsx)

if (!file.exists(paste0(sample_dir,"_quant.txt"))){
     quant <- featureCounts(files = print(paste0(MAP_DIR,"/",sample_id,"_marked_duplicates_STAR.bam")), annot.ext = repeat_gtf,isGTFAnnotationFile = T, isPairedEnd = TRUE, GTF.attrType = "gene_id",countMultiMappingReads = TRUE,primaryOnly = FALSE, fraction=TRUE,strandSpecific = strandness_fc)
     #colnames(quant)[i] <- sample_id
     write.table(quant$counts,print(paste0(quant_dir,"/",sample_id,"_quant.txt")),sep="\t")
}