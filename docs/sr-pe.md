# PERREO SR-PE

PERREO SR-PE is used for paired-end short-read RNA-seq data.

## Required Inputs

- Paired FASTQ files named as `<sample>_1.fastq` and `<sample>_2.fastq` inside each sample folder.
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
-adapt_r1
-adapt_r2
-trimming
-trimming_quality_threshold
-min_length_trim
-max_length_trim
-polya
-initial_trim_read1
-initial_trim_read2
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

SR-PE supports simple trimming with cutadapt and extra trimming for libraries where the kit adds GC nucleotides. Extra trimming uses `trimGC.py` after cutadapt and is only available for stranded paired-end data.

```bash
cutadapt -j "$threads" --pair-filter any -q "$trimming_quality","$trimming_quality" \
  -a "$adapt_r1" -A "$adapt_r2" \
  --trim-n -m "$min_length" -u "$initial_trim_read1" -U "$initial_trim_read2" "--$polya" \
  -o "${TRIM_DIR}/${sample_id}_trimmed_1.fastq" \
  -p "${TRIM_DIR}/${sample_id}_trimmed_2.fastq" \
  "$IN1" "$IN2" > cutadapt.log 2>&1
```

## Alignment

SR-PE uses STAR with multimapping enabled and a default mismatch threshold of `0.05`.

```bash
STAR --runThreadN "$threads" \
  --genomeDir "$GENOME_DIR" \
  --readFilesIn "$trimmed1" "$trimmed2" \
  --outSAMtype BAM SortedByCoordinate \
  --outFilterMultimapNmax 500 \
  --winAnchorMultimapNmax 500 \
  --outFilterMismatchNoverLmax "$mismatch_align"
```

## Quantification

Repeat RNA quantification is performed with `featureCounts` in paired-end mode, allowing multimapping reads and fractional assignment.
