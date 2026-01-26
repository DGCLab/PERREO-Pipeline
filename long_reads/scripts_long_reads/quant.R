args <- commandArgs(trailingOnly=TRUE)

sample_id <- args[[1]]
repeat_gtf <- args[[2]]
threads <- as.integer(args[[3]])
strandness <- args[[4]]
sample_dir <- args[[5]]

print(sample_id)
print(repeat_gtf)
print(threads)
print(strandness)
print(sample_dir)

if (strandedness=="unstranded"){
    strandness_fc = 0
}

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
     quant <- featureCounts(files = print(paste0(sample_dir,"/",sample_id,".bam")), annot.ext = repeat_gtf,isGTFAnnotationFile = T, isLongRead=T, isPairedEnd = FALSE, GTF.attrType = "gene_id",countMultiMappingReads = TRUE,primaryOnly = FALSE, fraction=TRUE,strandSpecific = strandness_fc, nthreads=threads)
     #colnames(quant)[i] <- sample_id
     write.table(quant$counts,print(paste0(sample_dir,"/",sample_id,"_quant.txt")),sep="\t")
     write.table(quant$stat, print(paste0(quant_dir,"/",sample_id,"_quant_stats.txt")),sep="\t")
}
