# PERREO

<p align="center">
  <img width="600" height="350" src="https://github.com/user-attachments/assets/f029a40a-dfc7-447c-b032-fe07a61f8a48">
</p>

<br>

<p align="center">
  <img src="https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white">
  <img src="https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white">
  <img src="https://img.shields.io/badge/Linux-FCC624?style=for-the-badge&logo=linux&logoColor=black">
  <img src="https://img.shields.io/badge/Conda-44A833?style=for-the-badge&logo=anaconda&logoColor=white">
</p>

# Software Requirements

To perform the entire analysis, it is necessary to create a conda environment with all required packages and programs. Users can achieve this by running the provided YAML file, which will install the complete environment automatically.

```bash
conda env create -f perreo.yml
```

# Data Preparation

The appropriate folder structure for the correct execution of the software is as follows:

```bash
PROJECT_FOLDER/
├─ samples
│ ├─ sample1
│ └─ sample2
├─ reference genome
├─ genome GTF
├─ repeat elements GTF
├─ samplesheet txt
└─ scripts
```

The "scripts" folder from this GitHub repository must be downloaded and placed inside the project folder.

The samplesheet structure should be as follows, regardless of the mode run:

| sample          | strandedness | condition |
|-----------------|--------------|-----------|
| SRR14506659     | reverse      | ESO       |
| SRR14506660     | reverse      | ESO       |
| SRR14506859     | reverse      | HC        |

If the user is unsure about the library strandedness, they can include a random type (forward/reverse) and decide later by applying specific software to infer strandedness. 
If data belong to more than two experimental conditions and the batch effect cause is known, an additional column called "batch" must be included in the samplesheet to indicate the origin batch.

For example: If samples come from different hospitals, this column should be included, and the hospital of origin must be indicated for each sample.

| sample          | strandedness | condition | batch |
|-----------------|--------------|-----------|-------|
| SRR14506659     | reverse      | ESO       | Hosp_A |
| SRR14506660     | reverse      | ESO       | Hosp_B |
| SRR14506859     | reverse      | HC        | Hosp_C |

Regarding the reference genome, users can choose any reference genome and annotations, as the software does not provide these files. Consequently, this pipeline is applicable to any organism with a sequenced and annotated genome.

For human experiments, although reference and annotations can be freely chosen, we recommend using the T2T genome and its corresponding annotation, as it is the most complete human reference with better-described repetitive regions.

Here, we provide a link to the T2T FASTA file and GTF annotations for use in this workflow: https://github.com/marbl/CHM13.

Another important point to consider is that both annotation files must indicate chromosomes consistently. For example, if genome annotations use "chrX" to indicate chromosome X, and repeat annotations use "X", the user must standardize one of them, as downstream analysis requires a uniform structure at certain points.

If the user needs to download other FASTA files and specific GTF annotations for T2T, the UCSC Table Browser platform can also be used: https://genome.ucsc.edu/cgi-bin/hgTables.

# Required Preanalysis

Before running the pipeline, perform FastQC and MultiQC to gather sufficient information for decisions regarding quality metrics before trimming reads and selecting whether to remove duplicates with MarkDuplicates.

To run MultiQC on FastQC outputs, use these command lines:

```bash
mkdir -p QC/fastqc_raw

# Process .fastq, .fq, .fastq.gz, and .fq.gz files
find . -type f \( -iname "*.fastq" -o -iname "*.fq" -o -iname "*.fastq.gz" -o -iname "*.fq.gz" \) -print0 \
| xargs -0 -n 1 -P 4 fastqc -t 2 -o QC/fastqc_raw
```

# Running PERREO via User-Friendly Interface

To run PERREO from a graphical interface, download the "app.R" script and place it in the general directory where the pipeline will be executed. Then, run the following code:

```bash
export PATH="/home/user/project_directory:$PATH"
Rscript app.R
```

Once the script is executed, the interface will open in the web browser, allowing users to select parameters and upload required files in a user-friendly manner. During execution, the program will print messages related to the process in the terminal, where users can monitor progress.



# Workflow Summary

The flow diagram illustrates the different steps included in this pipeline.

![PERREO](https://github.com/user-attachments/assets/3385d410-6be1-4de8-97a1-d07a007d40af)

# PERREO Modes

PERREO offers three modes depending on the sequencing technology used:

- PERREO SR-SE: For RNA-seq data generated with single-end short-read technology.
- PERREO SR-PE: For RNA-seq data generated with paired-end short-read technology.
- PERREO LR: For direct RNA-seq data generated with Nanopore long-read technology.

<img width="2437" height="875" alt="logos" src="https://github.com/user-attachments/assets/5240009a-2645-4fd1-9fff-7e5ced4cf3b5" />


# General Requirements

These arguments are required for the correct execution of the corresponding modes:

| Input/Argument          | PERREO SR-SE | PERREO SR-PE | PERREO LR |
|-------------------------|--------------|--------------|-----------|
| Single-end FASTQ files  | X            |              | X         |
| Paired-end FASTQ files  |              | X            |           |
| Reference genome        | X            | X            | X         |
| Genome GTF              | X            | X            |           |
| Repeats GTF             | X            | X            | X         |
| Batch effect            | X            | X            | X         |
| Method                  | X            | X            | X         |
| Remove duplicates       | X            | X            |           |

Additional arguments to consider, though not mandatory as they have default values:

| Argument                          | Default |
|-----------------------------------|---------|
| Threads                           | 8       |
| RAM                               | 32      |
| Trimming                          | simple (not needed for PERREO SR-SE & LR) |
| Mismatch align                    | 0.05    |
| Trimming quality threshold        | 30      |
| Minimum read length (for trimming)| 16      |
| Maximum read length (for trimming)|         |
| Initial nucleotides trimming R1 (paired-end) | 0 |
| Initial nucleotides trimming R2 (paired-end) | 0 |
| Initial nucleotides trimming (single-end) | 0 |
| K_num                             | 2       |
| log2FC                            | 1       |
| FDR                               | 0.05    |
| Prediction model                  | no      |
| Batch                             | no      |

For the LR mode, trimming parameters such as trimming, trimming_quality_threshold, min_length_trim, and max_length_trim are not considered, as trimming for long reads should be performed during basecalling before obtaining FASTQ files.

In this documentation, we first describe the trimming and alignment processes, which are specific to each PERREO mode. Then, we detail the downstream analysis, which is common across the three modes, except for the quantification step, where specific featureCounts arguments must be used for each case.


# PERREO SR-PE and SR-SE

SR-SE and SR-PE analyses are very similar, with differences only in the trimming and alignment steps. The rest of the analysis is common between these two modes.

The arguments for SR-PE and SR-SE analyses are nearly identical. However, there are specific differences, which we will describe separately for each PERREO mode.

For SR-PE:

```text
-sample_list                  Sample sheet with sample, strandedness, condition, and batch (if necessary).
-reference_genome             Genome file in FASTA format.
-genome_gtf                   Genome annotations in GTF format.
-repeat_gtf                   Repeat annotations in GTF format. It must contain a "repeat_class" column to study the type of repetitive elements identified in the analysis.
-threads                      Number of threads used for the process (default: 8).
-ram                          RAM memory used for MarkDuplicates in GB (default: 32).
-adapt_r1                     Adapter 1 sequence. If reads are already trimmed, ignore this argument and set to "" for high-quality read trimming.
-adapt_r2                     Adapter 2 sequence. If reads are already trimmed, ignore this argument and set to "" for high-quality read trimming.
-trimming                     simple/extra. Select "extra" if the kit adds extra GC nucleotides (default: simple).
-trimming_quality_threshold   Minimum quality permitted for reads to be kept after trimming (default: 30).
-min_length_trim              Minimum read length to retain after trimming (default: 16).
-max_length_trim              Maximum read length permitted for reads to be kept after trimming (default: none).
-polya                        If polyA tails need removal (default: not applied).
-initial_trim_read1           Number of initial nucleotides to trim in R1 files after adapter trimming (default: 0).
-initial_trim_read2           Number of initial nucleotides to trim in R2 files after adapter trimming (default: 0).
-mismatch_align               Percentage of mismatches permitted for reads during alignment (default: 0.05).
-project_name
-remove_duplicates            yes/no. Generally not recommended unless duplication proportion is very high.
-batch                        yes/no. If yes, and batch column is included in the sample sheet, batch effect will be reduced based on that column. If yes but no batch column, RUVg will reduce undesired variability.
-method                       edgeR/DESeq2.
-k_num                        K number for RUVg-based batch effect reduction. Increase with dataset size; 2-3 factors usually suffice (default: 2).
-log2FC                       Log2-transformed fold-change threshold for Differential Expression Analysis (default: 1).
-FDR                          Adjusted p-value threshold for Differential Expression Analysis (default: 0.05).
-prediction_model             yes/no (default: no). When activated, designs prediction models only if samples > 40.
-positive_class               The experimental condition to be the positive class in the prediction model.
```

For SR-SE:

```text
-sample_list                  Sample sheet with sample, strandedness, condition, and batch (if necessary).
-reference_genome             Genome file in FASTA format.
-genome_gtf                   Genome annotations in GTF format.
-repeat_gtf                   Repeat annotations in GTF format. It must contain a "repeat_class" column to study the type of repetitive elements identified in the analysis.
-threads                      Number of threads used for the process (default: 8).
-ram                          RAM memory used for MarkDuplicates in GB (default: 32).
-adapter                      Adapter sequence. If reads are already trimmed, ignore this argument and set to "" for high-quality read trimming.
-trimming_quality_threshold   Minimum quality permitted for reads to be kept after trimming (default: 30).
-min_length_trim              Minimum read length to retain after trimming (default: 16).
-max_length_trim              Maximum read length permitted for reads to be kept after trimming (default: none).
-polya                        If polyA tails need removal (default: not applied).
-initial_trim_read            Number of initial nucleotides to trim in FASTQ file after adapter trimming (default: 0).
-mismatch_align               Percentage of mismatches permitted for reads during alignment (default: 0.05).
-project_name
-remove_duplicates            yes/no. Generally not recommended unless duplication proportion is very high.
-batch                        yes/no. If yes, and batch column is included in the sample sheet, batch effect will be reduced based on that column. If yes but no batch column, RUVg will reduce undesired variability.
-method                       edgeR/DESeq2.
-k_num                        K number for RUVg-based batch effect reduction. Increase with dataset size; 2-3 factors usually suffice (default: 2).
-log2FC                       Log2-transformed fold-change threshold for Differential Expression Analysis (default: 1).
-FDR                          Adjusted p-value threshold for Differential Expression Analysis (default: 0.05).
-prediction_model             yes/no (default: no). When activated, designs prediction models only if samples > 40.
-positive_class               The experimental condition the user wants to be the positive class in the prediction model.
```


## Trimming

In this step, there are two main options: simple trimming with cutadapt, and a more complex trimming performed first with cutadapt and then remove additional GC nucleotides added by specific sequencing kits.

If adapters have been removed previously and trimming is only needed to remove low-quality reads, include the -adapt_r1, -adapt_r2, and -adapter arguments as indicated in the previous sections while running the pipeline.

To remove polyA tails, include the -polya argument as follows: -polya polya.


Single-end cutadapt trimming:

```bash
cutadapt -j "$threads" -q "$trimming_quality","$trimming_quality" \
        -a "$adapter" --trim-n -m "$min_length" -u "$initial_trim_read" \
        -o "${TRIM_DIR}/${sample_id}_trimmed.fastq" "--$polya" \
        "$IN" > cutadapt.log 2>&1
```

Paired-end cutadapt trimming:

```bash
cutadapt -j "$threads" --pair-filter any -q "$trimming_quality","$trimming_quality" \
        -a "$adapt_r1" -A "$adapt_r2" \
        --trim-n -m "$min_length" -u "$initial_trim_read1" -U "$initial_trim_read2" "--$polya" \
        -o "${TRIM_DIR}/${sample_id}_trimmed_1.fastq" \
        -p "${TRIM_DIR}/${sample_id}_trimmed_2.fastq" \
        "$IN1" "$IN2" > cutadapt.log 2>&1
```

## Alignment

STAR aligner parameters are set by default to allow multimapping and remove reads with more than 5% mismatches. This parameter can be adjusted based on the user's experimental design and conditions.

```bash
STAR --runThreadN $threads \
      --genomeDir "$GENOME_DIR" \
      --readFilesIn "$trimmed1" "$trimmed2" \ # for single-end data, only "$trimmed"
      --outFileNamePrefix "$MAP_DIR/${sample_id}_" \
      --outSAMtype BAM SortedByCoordinate \
      --outFilterMultimapNmax 500 \
      --winAnchorMultimapNmax 500 \
      --outFilterMismatchNoverLmax 0.05 \
      --outSAMmultNmax 500 \
      --outMultimapperOrder Random \
      --runRNGseed 42 \
      --outSAMattributes NH HI AS nM NM MD
```

## Duplicates Analysis

Duplicate removal is generally not recommended in RNA-seq data analysis. However, in cases where the duplication percentage is excessively high, removal may be an option.

## Quantification

FeatureCounts performs feature quantification, allowing multimapping read counts and fractions. It uses the strandedness indicated in the sample sheet.

A quantification statistics file is exported for each sample, detailing assigned and unassigned reads.

The following code is run for data from single-end and paired-end short-read sequencing technologies:

Single-end:

```R
quant <- featureCounts(files = print(paste0(MAP_DIR,"/",sample_id,"_marked_duplicates_STAR.bam")), annot.ext = repeat_gtf,
         isGTFAnnotationFile = T, isPairedEnd = FALSE, GTF.attrType = "gene_id", countMultiMappingReads = TRUE, primaryOnly = FALSE,
         fraction = TRUE, strandSpecific = strandness_fc)
```

Paired-end:

```R
quant <- featureCounts(files = print(paste0(MAP_DIR,"/",sample_id,"_marked_duplicates_STAR.bam")), annot.ext = repeat_gtf,
         isGTFAnnotationFile = T, isPairedEnd = TRUE, GTF.attrType = "gene_id", countMultiMappingReads = TRUE, primaryOnly = FALSE,
         fraction = TRUE, strandSpecific = strandness_fc)
```


# PERREO LR

The required arguments for this mode are as follows:

```text
-sample_list           Sample sheet with sample, strandedness, condition, and batch (if necessary).
-reference_genome      Genome file in FASTA format.
-repeat_gtf            Repeat annotations in GTF format. It must contain a "repeat_class" column to study the type of repetitive elements identified in the analysis.
-threads               Number of threads used for the process (default: 8).
-project_name
-batch_effect          yes/no. If yes, and batch column is included in the sample sheet, batch effect will be reduced based on that column. If yes but no batch column, RUVg will reduce undesired variability.
-method                edgeR/DESeq2.
-k_num                 K number for RUVg-based batch effect reduction. Increase with dataset size; 2-3 factors usually suffice (default: 2).
-log2FC                Log2-transformed fold-change threshold for Differential Expression Analysis (default: 1).
-FDR                   Adjusted p-value threshold for Differential Expression Analysis (default: 0.05).
-prediction_model      yes/no (default: no). When activated, designs prediction models only if samples > 40.
-positive_class        The experimental condition to be the positive class in the prediction model.
```

In this case, trimming parameters are not required, as adapter and barcode removal is performed during basecalling, and quality control is conducted using NanoPlot.

## Alignment

Long reads are aligned against the reference genome using minimap2 with -ax splice -uf and -k14 parameters. Multimapping is allowed with -N 100. It is important to set -ax splice for splice-aware alignment and -uf to assume forward orientation relative to the transcribed strand.

```bash
minimap2 -t 14 -ax splice -uf -k14 -p 0.8 -N 100 "$CWD/genome_index.mmi" "$SAMPLE_DIR/${sample_id}.fastq" > "$SAMPLE_DIR/${sample_id}.sam"
```

## Quantification

For data from long-read technology, the code is nearly identical, with the long-read argument activated:

```R
quant <- featureCounts(files = print(paste0(sample_dir,"/",sample_id,".bam")), annot.ext = repeat_gtf, isGTFAnnotationFile = T,
         isLongRead = T, isPairedEnd = FALSE, GTF.attrType = "gene_id", countMultiMappingReads = TRUE, primaryOnly = FALSE, fraction = TRUE,
         strandSpecific = strandness_fc)
```

# Downstream Analysis

Differential expression analysis, coexpression analysis, transcriptome assembly, and prediction model design are common across the three PERREO modes.

## Differential Expression Analysis

Statistical analysis is performed using DESeq2 or edgeR with default thresholds: abs(log2FC) > 1 and FDR < 0.05. DESeq2 is robust and conservative, while edgeR is sensitive for low-expressed genes due to flexible dispersion estimation. If batch effect = yes, RUVg is applied if no batch column exists in the sample sheet. If present, it is included in the formula. Various plots and files are exported:

1. Repeat RNA counts represented with violin plots.
2. PCA (if batch effect is corrected, PCA before and after correction are exported).
3. Fold-change bar plots of differentially expressed features (DEFs) for each contrast.
4. Classification of all identified repeat RNAs by repeat classes.
5. Classification of differentially expressed repeat RNAs by repeat classes.
6. Expression matrix in CSV format.
7. Heatmap of DEFs in each contrast.
8. Volcano plot for each contrast.

## Coexpression Analysis

The WGCNA R package generates coexpression networks from the expression matrix. Power is automatically assigned when R = 0.9, and the script selects the three modules with the highest correlation to any experimental condition. It exports:

1. Information about nodes and edges of the three modules for Cytoscape networks.
2. Correlation heatmap between modules and experimental conditions.
3. Correlation heatmap between modules and samples.
4. Repeat RNAs and their modules in TXT format.
5. Modules dendrogram.

## Transcriptome Assembly

Stringtie2 generates a GTF file for transcriptome assembly of each sample using a combined GTF of genomic and repeat annotations as reference. Strandedness must be considered.

```bash
stringtie "$BAM" \
        -L \
        -p "$threads" \
        -G "$combined_annotation" \
        -o "$OUTPUT_GTF" \
        --fr \ # rf if strandedness=reverse
```

Then, the pipeline performs a preliminary analysis mapping each sample's annotation to GTFs containing exonic regions and repetitive elements to identify hybrid transcripts composed of exonic and repeat fragments.

## Prediction Models

GLMnet and Random Forest algorithms generate prediction models to assess the classifying potential of repeat RNAs. It is important to define the positive class and provide sufficient samples.

A common error is:

```text
Error in roc.default(bin_truth, prob_df[[cls]], levels = rev(levels(bin_truth)),  :
  No case observation.
Calls: eval_and_export ... mutate -> ovr_roc_curves -> <Anonymous> -> roc.default
Execution halted
```

This may arise from insufficient samples per class or severe imbalance. The pipeline applies balancing before cross-validation to ensure representation. If persistent, verify sufficient samples per class relative to folds, and consider merging or excluding underrepresented classes. In such cases, remove specific samples from the samplesheet and rerun if prior steps are complete. This way, samples are included in DEA, coexpression, and assembly but not prediction models.

## PDF Report

The generated report includes the plots mentioned in the Differential Expression section.

<img width="6000" height="4500" alt="VolcanoPlot" src="https://github.com/user-attachments/assets/fbcd3692-6819-4f7e-a95f-46432376ff58" />

<img width="2400" height="1800" alt="repetitive_counts_violin_box" src="https://github.com/user-attachments/assets/c7bde717-9a58-43da-89f2-2b35d2c12b94" />

<img width="2400" height="1800" alt="pca" src="https://github.com/user-attachments/assets/a16a83d9-48bf-4048-91bc-e6553f1f734f" />

<img width="6000" height="4500" alt="Classification_All" src="https://github.com/user-attachments/assets/37aa337a-d45d-4778-be2f-90bb57105c25" />

<img width="6000" height="4500" alt="BarPlotUpDown" src="https://github.com/user-attachments/assets/69d56e63-f1c4-4996-b6a4-8ff4072e35c1" />



