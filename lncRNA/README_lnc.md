# LncRNAs analysis<br>
Many long non-coding RNAs (lncRNAs) are also active in the cell transcriptome. They are not included in repeatmasker annotations tracks, but they can be found in GENCODE specific annotations.<br>

Here we propose a comprehensive pipeline to study the lncRNAs expression in different types of samples. The pipeline is composed of the same steps as PERREO, but changing some specific arguments specialized in detecting long non-coding features. The workflow is designed to run the bioinformatic analysis with data derived from paired-end sequencing approaches. <br>

In contrast to PERREO run code line, the user has to provide a long non-coding elements annotation file (from GENCODE, for example) instead of a repeatmasker track annotations GTF. Additionally, the pipeline run STAR and featureCounts with different parameters.<br>

## STAR alignment<br>
For reads mapping, parameters related to multimappers management are dramatically reduced as lncRNAs must not map to multiple sites as repetitive elements do. The rest of the code line is mantained.

```bash
STAR --runThreadN $threads \
      --genomeDir "$GENOME_DIR" \
      --readFilesIn "$trimmed1" "$trimmed2" \
      --outFileNamePrefix "$MAP_DIR/${sample_id}_" \
      --outSAMtype BAM SortedByCoordinate \
      --outFilterMultimapNmax 10 \
      --winAnchorMultimapNmax 50 \
      --outFilterMismatchNoverLmax $mismatch_align \
      --outSAMmultNmax 10 \
      --outMultimapperOrder Random \
      --runRNGseed 42 \
      --outSAMattributes NH HI AS nM NM MD
```

## featureCounts quantification<br>
Regarding quantification, multimapping reads are not taken into account, as well as fractioning is also deactivated in order to only detect lncRNAs that truly map to specific genomic regions. 
```bash
quant <- featureCounts(files = print(paste0(MAP_DIR,"/",sample_id,"_marked_duplicates_STAR.bam")), annot.ext = lnc_gtf,isGTFAnnotationFile = T, isPairedEnd = TRUE, GTF.attrType = "gene_id",countMultiMappingReads = FALSE,primaryOnly = FALSE, fraction=FALSE,strandSpecific = strandness_fc, nthreads=threads)

```


