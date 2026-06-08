# PERREO SR-SE

PERREO SR-SE is used for single-end short-read RNA-seq data.

## Required Inputs

- FASTQ files named as `<sample>.fastq` inside each sample folder.
- Samplesheet with `sample`, `strandedness`, and `condition`.
- Reference genome FASTA.
- Genome GTF.
- Repeat GTF with a `repeat_class` attribute.

## Main Arguments

```text
-sample_list
-reference_genome
-genome_gtf
-repeat_gtf
-threads
-ram
-adapter
-trimming_quality_threshold
-min_length_trim
-max_length_trim
-polya
-initial_trim_read
-mismatch_align
-project_name
-remove_duplicates
-batch
-method
-k_num
-log2FC
-FDR
-prediction_model
-positive_class
```

## Trimming

SR-SE uses cutadapt for adapter removal, quality filtering, optional initial base trimming, and optional polyA removal.

```bash
cutadapt -j "$threads" -q "$trimming_quality","$trimming_quality" \
  -a "$adapter" --trim-n -m "$min_length" -u "$initial_trim_read" \
  -o "${TRIM_DIR}/${sample_id}_trimmed.fastq" "--$polya" \
  "$IN" > cutadapt.log 2>&1
```

## Alignment

SR-SE uses STAR with the same multimapping and mismatch settings as SR-PE.

```bash
STAR --runThreadN "$threads" \
  --genomeDir "$GENOME_DIR" \
  --readFilesIn "$trimmed" \
  --outSAMtype BAM SortedByCoordinate \
  --outFilterMultimapNmax 500 \
  --winAnchorMultimapNmax 500 \
  --outFilterMismatchNoverLmax "$mismatch_align"
```

## Quantification

Repeat RNA quantification is performed with `featureCounts` in single-end mode, allowing multimapping reads and fractional assignment.
