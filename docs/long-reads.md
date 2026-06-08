# PERREO LR

PERREO LR is used for Nanopore direct RNA-seq long-read data.

## Required Inputs

- FASTQ files named as `<sample>.fastq` inside each sample folder.
- Samplesheet with `sample`, `strandedness`, and `condition`.
- Reference genome FASTA.
- Repeat GTF with a `repeat_class` attribute.

Long-read mode does not require genome GTF for quantification, but transcriptome assembly steps use combined annotations when available in the workflow.

## Main Arguments

```text
-sample_list
-reference_genome
-repeat_gtf
-threads
-project_name
-batch_effect
-method
-k_num
-log2FC
-FDR
-prediction_model
-positive_class
```

## Quality Control

Adapter and barcode trimming should be performed during basecalling before PERREO receives FASTQ files. Long-read quality control is performed with NanoPlot.

## Alignment

Long reads are aligned with minimap2 using splice-aware settings:

```bash
minimap2 -t 14 -ax splice -uf -k14 -p 0.8 -N 100 \
  "$CWD/genome_index.mmi" "$SAMPLE_DIR/${sample_id}.fastq" \
  > "$SAMPLE_DIR/${sample_id}.sam"
```

The BAM file is then sorted and indexed with samtools.

## Quantification

Repeat RNA quantification is performed with `featureCounts` using long-read mode:

```r
featureCounts(
  files = paste0(sample_dir, "/", sample_id, ".bam"),
  annot.ext = repeat_gtf,
  isGTFAnnotationFile = TRUE,
  isLongRead = TRUE,
  isPairedEnd = FALSE,
  GTF.attrType = "gene_id",
  countMultiMappingReads = TRUE,
  primaryOnly = FALSE,
  fraction = TRUE,
  strandSpecific = strandness_fc
)
```
