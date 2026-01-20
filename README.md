# PERREO<br> 

<img width="600" height="350" alt="logo_github" src="https://github.com/user-attachments/assets/f029a40a-dfc7-447c-b032-fe07a61f8a48" />


# Software requirements:<br> 

It is necessary to create a conda environment with all the required packages and programs to perform the whole analysis. To do it, the user can run this yaml file. In this way, the complete environment will be installed automatically. 
```bash
conda env create -f perreo.yml

```

# Data preparation<br> 

The suitable folders structure for the correct performance of the software should be the following:<br> 

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

The folder "scripts" located in this github must be downloaded and included inside the project folder.<br> 


The samplesheet structure should be like this regardless of the mode run:<br>


| sample          | strandedness | condition |
|-----------------|--------------|-----------|
| SRR14506659     | reverse      | ESO       |
| SRR14506660     | reverse      | ESO       |
| SRR14506859     | reverse      | HC        |

If the user does not know the library strandedness, he/she can include a random type (forward/reverse) and decide after alignment applying specific software to infer the strandedness, like infer_experiment.py or check_for_strandedness. However, in this case, the user cannot perform extra trimming as trimGC.py script requires the strandedness to trim the reads. <br>

In case data belong to more than two experimental conditions and batch effect cause is known, another column called "batch" must be included in the samplesheet indicating the origin batch. <br>
<br>
For example: If samples come from different hospitals, this column should be included and the hospital of origin has to be indicated for each sample.<br>

| sample          | strandedness | condition | batch |
|-----------------|--------------|-----------|-------|
| SRR14506659     | reverse      | ESO       |Hosp_A |
| SRR14506660     | reverse      | ESO       |Hosp_B |
| SRR14506859     | reverse      | HC        |Hosp_C |

<br>
Regarding the reference genome, the user can decide which reference genome and which annotations to use as the software does not provide any of these files. Consequently, this pipeline is applicable to any organism whose genome is sequenced and annotated.<br> 
<br>
For human experiments, although reference and annotations can be freely chosen, we recommend to use T2T genome and its corresponding annotation as it is the most complete human reference, where repetitive regions are better described.<br>
<br>
Here we provide a link to the T2T fasta file and GTF annotations to use in this workflow: https://github.com/marbl/CHM13.<br>
<br>
Another important point to take into account is that both annotations files must indicate chromosomes in the same way. For example, if genome annotations use the structure "chrX" to indicate chromosome X, and repeat annotationes refers to it using "X", the user has to change one of them because downstream analysis requires only one possible structure at certain point.<br>
<br>
In case the user needs to download other fasta files and specific GTF annotations of T2T, the UCSC table browser platform can also be used: https://genome.ucsc.edu/cgi-bin/hgTables.<br>

# Required preanalysis<br>

Before running the pipeline you should perform fastqc and multiqc to have the enough information to make some decisions with respect to quality metrics before trimming the reads and selecting or not the remove duplicates option of MarkDuplicates.<br>


To run multiqc on fastqc outputs, the user can use these command lines:

```bash
mkdir -p QC/fastqc_raw

# Process .fastq, .fq, .fastq.gz y .fq.gz
find . -type f \( -iname "*.fastq" -o -iname "*.fq" -o -iname "*.fastq.gz" -o -iname "*.fq.gz" \) -print0 \
| xargs -0 -n 1 -P 4 fastqc -t 2 -o QC/fastqc_raw
```

# Running PERREO on user-friendly interface<br>
To run PERREO from a graphical interface the user has to download the script "app.R" and store it in the general directory where the pipeline is going to be run. Then, the user has to run the following code:

```bash
export PATH="/home/user/project_directory:$PATH"
Rscript app.R
```
Once the app.R script is run, the interface will be open in the internet explorer and the user will be able to select all the parameters and upload the required files in a user-friendly manner. When the workflow is executed, the program will print different messages related to the execution in the terminal, where the user will be able to check the process.



# Workflow summary<br>

The flow diagram describes the different steps taken into account in this pipeline.<br>
<br>
![PERREO](https://github.com/user-attachments/assets/3385d410-6be1-4de8-97a1-d07a007d40af)

<br>

# PERREO modes<br>

PERREO provides three models that can be run depending on the sequencing technology used:<br>

PERREO SR-SE: For RNA-seq data generated with single-end short reads technology.<br>
PERREO SR-PE: For RNA-seq data generated with paired-end short-reads technology.<br>
PERREO LR: For direct RNA-seq data generated with Nanopore long-reads technology.<br>

<img width="2437" height="875" alt="logos" src="https://github.com/user-attachments/assets/5240009a-2645-4fd1-9fff-7e5ced4cf3b5" />


# General requirements<br>

These arguments are required for the correct performance of the corresponding modes:

| Input/argument          | PERREO SR-SE | PERREO SR-PE | PERREO LR |
|-----------------|--------------|-----------|-------|
| Single-end fastq files     | X      |        |X |
| Paired-end fastq files     |       | X       | |
| Reference genome     | X      | X        |X |
| Genome GTF     | X      | X       |  |
| Repeats GTF     | X      | X       |X |
| Batch effect     | X      | X       |X |
| Method     | X      | X       |X |
| Remove duplicates     | X      | X       | |


Then, there are some arguments that have to be taken into account although it is not mandatory to mention them in the command line because they have default values:
| Arguments          | Default |
|-----------------|--------------|
| Threads     | 8      |
| Trimming      |  simple (not needed for PERREO LR)     | 
| Mismatch align     | 0.05     | 
| Trimming quality threshold     | 30     | 
| Minimum reads length (for trimming)     | 16     | 
| Maximum reads length (for trimming)     |      | 
| Initial nucleotides trimming R1 (paired-end)     |   0   | 
| Initial nucleotides trimming R2 (paired-end)    |   0   | 
| Initial nucleotides trimming (single-end)     |   0   | 
| K_num     | 2     | 
| log2FC     | 1     | 
| FDR     | 0.05     | 
| Prediction model     | no     |
| batch     | no     | 



For SR-LR mode, the parameters trimming_type, trimming_quality_threshold, min_length_trim and max_length_trim are not considered as trimming for long-reads should be performed during the basecalling before obtaining the fastq files.<br>




In the this documentation, we firstly describe the trimming and alignment process, which is specific for each PERREO mode. Then, we describe in detail how the downstream analysis is performed, which is common between the three modes, excepting the quantification step where specific featureCounts arguments must be included in each case.


# PERREO SR-PE and SR-SE
SR-SE and SR-PE analyses are very similar as the only differences are located in the trimming and alignment steps. The rest of the analysis is common  between these two modes.<br> 
The arguments for SR-PE and SR-SE analysis are practically the same. However, there exists specific differences and we will show them separately by PERREO mode. <br> 

For SR-PE:<br> 
```text
-sample_list                  Sample sheet with sample, strandedness, condition and batch (if necessary).
-reference_genome             Genome file in fasta.
-genome_gtf                   Genome annotations in GTF format.
-repeat_gtf                   Repeat annotations in GTF format. It must contain a "repeat_class" column in order to study the type
                              of repetitive elements identified in the analysis.
-threads                      Number of threads used for the process (default: 8)
-adapt_r1                     Adapter 1 sequence. If reads are already trimmed, the user can ignore this argument and trimming will
                              be also performed to obtain high quality reads by indicating "-adapt_r1 """.
-adapt_r2                     Adapter 2 sequence. If reads are already trimmed, the user can ignore this argument and trimming will
                              be also performed to obtain high quality reads by indicating "-adapt_r2 """.
-trimming                     simple/extra. The second must be selected if it is know that the kit used adds extra
                              GC nucleotides (default: trimming_simple).
-trimming_quality_threshold   Minimum quality permitted for reads to be kept after trimming (default: 30).
-min_length_trim              Minimum reads length to not discard them after trimming (default: 16).
-max_length_trim              Maximum reads length permitted for reads to be kept after trimming (default: ).
-polya                        If polyA tails have to be removed (default: not applied).
-initial_trim_read1           Number of initial nucleotides that must be trimmed in R1 files after specific adapter trimming (default: 0).
-initial_trim_read2           Number of initial nucleotides that must be trimmed in R2 files after specific adapter trimming (default: 0).
-mismatch_align               Percentage of mismatches permitted for reads to be kept during alignment (default: 0.05).
-project_name
-remove_duplicates            yes/no. Generally it is not recommended to discard duplicates unless the proportion is really high.
-batch                        yes/no. If yes, and batch column is included in the sample sheet, the batch effect will be reduced
                              based on that column. If yes but batch column does not exist, RUVg will reduce the undesired variability.
-method                       edgeR/DESeq2. 
-k_num                        K number for RUVg-based batch effect reduction. As far as the dataset size increases, it is recommended
                              to also increase the k number (default: 2).
-log2FC                       Log2-transformed fold-change threshold for Differential Expression Analysis (default: 1).
-FDR                          Adjusted p-value threshold for Differential Expression Analysis (default: 0.05).
-prediction_model             yes/no (default: no). When activated, it only will design prediction models in case the number of samples
                              is higher than 40.
-positive_class               The experimental condition the user wants to be the positive class in the prediction model.
```
For SR-SE:<br> 

```text
-sample_list                  Sample sheet with sample, strandedness, condition and batch (if necessary).
-reference_genome             Genome file in fasta.
-genome_gtf                   Genome annotations in GTF format.
-repeat_gtf                   Repeat annotations in GTF format. It must contain a "repeat_class" column in order to study the type
                              of repetitive elements identified in the analysis.
-threads                      Number of threads used for the process (default: 8)
-adapter                      Adapter sequence. If reads are already trimmed, the user can ignore this argument and trimming will
                              be also performed to obtain high quality reads by indicating "-adapter """.
-trimming                     simple/extra. The second must be selected if it is know that the kit used adds extra
                              GC nucleotides (default: trimming_simple).
-trimming_quality_threshold   Minimum quality permitted for reads to be kept after trimming (default: 30).
-min_length_trim              Minimum reads length to not discard them after trimming (default: 16).
-max_length_trim              Maximum reads length permitted for reads to be kept after trimming (default: ).
-polya                        If polyA tails have to be removed (default: not applied).
-initial_trim_read            Number of initial nucleotides that must be trimmed in fastq file after specific adapter trimming (default: 0).
-mismatch_align               Percentage of mismatches permitted for reads to be kept during alignment (default: 0.05).
-project_name
-remove_duplicates            yes/no. Generally it is not recommended to discard duplicates unless the proportion is really high.
-batch                        yes/no. If yes, and batch column is included in the sample sheet, the batch effect will be reduced
                              based on that column. If yes but batch column does not exist, RUVg will reduce the undesired variability.
-method                       edgeR/DESeq2. 
-k_num                        K number for RUVg-based batch effect reduction. As far as the dataset size increases, it is recommended
                              to also increase the k number (default: 2).
-log2FC                       Log2-transformed fold-change threshold for Differential Expression Analysis (default: 1).
-FDR                          Adjusted p-value threshold for Differential Expression Analysis (default: 0.05).
-prediction_model             yes/no (default: no). When activated, it only will design prediction models in case the number of samples
                              is higher than 40.
-positive_class               The experimental condition the user wants to be the positive class in the prediction model.

```


## Trimming
In this step there are two main options: simple trimming with cutadapt and a more complex trimming performed firstly with cutadapt and then with trimGC.py script in order to remove additional GC nucleotides added by specific sequencing kits. <br>
<br>
If adapters had been removed previously and the trimming step should only be taken into account to remove low quality reads, user must include -adapt_r1, -adapt_r2 and -adapter arguments as indicated in the previous boxes while running the pipeline.<br>
In case the user also want to remove polyA tails, he/she must include that -polya argument as following: -polya polya.<br>
If the library used was unstranded, even if the user selects extra trimming, the pipeline will perform simple trimming, as trimGC.py script only can be executed on data derived from forward and reverse libraries. In this way, simple trimming is automatically performed on unstranded library-based data.

Single-end cutadapt trimming:
```bash
 cutadapt -j "$threads" -q "$trimming_quality","$trimming_quality" \
        -a "$adapter" --trim-n -m "$min_length" -u "$initial_trim_read" \
        -o "${TRIM_DIR}/${sample_id}_trimmed.fastq" "--$polya" \
        "$IN"  > cutadapt.log 2>&1
```

Paired-end cutadapt trimming:
```bash
   cutadapt -j "$threads" --pair-filter any -q "$trimming_quality","$trimming_quality" \
        -a "$adapt_r1" -A "$adapt_r2" \ 
        --trim-n -m "$min_length" -u "$initial_trim_read1"  -U "$initial_trim_read2" "--$polya" \
        -o "${TRIM_DIR}/${sample_id}_trimmed_1.fastq" \ 
        -p "${TRIM_DIR}/${sample_id}_trimmed_2.fastq" \ 
        "$IN1" "$IN2" > cutadapt.log 2>&1
```

## Alignment
STAR aligner parameters are indicated by default allowing the multimapping and removing reads with more than a 5% of mismatches. This parameter can be changed depending on the user experimental design and conditions.<br>

```bash
STAR --runThreadN $threads \
      --genomeDir "$GENOME_DIR" \
      --readFilesIn "$trimmed1" "$trimmed2" \ #for single-end data only "$trimmed"
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

## Duplicates analysis
Duplicateds removal is not recommended generally in RNA-seq data analysis. However, there are situations where the percentage of duplications is too high and remove them is a real option.

## Quantification
FeatureCounts performs features quantification allowing multimapping reads counts and fraction. It uses the strandedness indicated in the sample sheet.<br>

The following code line is run for data obtained from single-end and paired-end short reads sequencing techonologies:<br>

Single-end:<br>

```R
quant <- featureCounts(files = print(paste0(MAP_DIR,"/",sample_id,"_marked_duplicates_STAR.bam")), annot.ext = repeat_gtf,
         isGTFAnnotationFile = T, isPairedEnd = FALSE, GTF.attrType = "gene_id",countMultiMappingReads = TRUE,primaryOnly = FALSE,
         fraction=TRUE,strandSpecific = strandness_fc)
```

Paired-end:<br>

```R
quant <- featureCounts(files = print(paste0(MAP_DIR,"/",sample_id,"_marked_duplicates_STAR.bam")), annot.ext = repeat_gtf,
         isGTFAnnotationFile = T, isPairedEnd = TRUE, GTF.attrType = "gene_id",countMultiMappingReads = TRUE,primaryOnly = FALSE,
         fraction=TRUE,strandSpecific = strandness_fc)
```


# PERREO LR <br> 
The required arguments for this mode are the following:<br>

```text
-sample_list           Sample sheet with sample, strandedness, condition and batch (if necessary).
-reference_genome      Genome file in fasta.
-repeat_gtf            Repeat annotations in GTF format. It must contain a "repeat_class" column in order to study the type
                       of repetitive elements identified in the analysis.
-threads               Number of threads used for the process (default: 8)
-project_name
-batch_effect          yes/no. If yes, and batch column is included in the sample sheet, the batch effect will be reduced
                       based on that column. If yes but batch column does not exist, RUVg will reduce the undesired variability.
-method                edgeR/DESeq2. 
-k_num                 K number for RUVg-based batch effect reduction. As far as the dataset size increases, it is recommended
                       to also increase the k number (default: 2).
-log2FC                Log2-transformed fold-change threshold for Differential Expression Analysis (default: 1).
-FDR                   Adjusted p-value threshold for Differential Expression Analysis (default: 0.05).
-prediction_model      yes/no (default: no). When activated, it only will design prediction models in case the number of samples
                       is higher than 40.
-positive_class        The experimental condition the user wants to be the positive class in the prediction model.

```

In this case, trimming parameters are not required as adapters and barcodes removal is carried out during the basecalling step and quality control is carried out using NanoPlot.<br>

## Alignment
Long reads are aligned against the reference genome using minimap2 with -ax splice -uf and -k14 parameters. Multimapping is also allowed indicating -N 100 parameter. It is important to set -ax splice, in order to perform splice-aware alignment, and -uf argument to try to align reads assuming it follows the forward orientation relative to the transcribed strand.

```bash
minimap2 -t 14 -ax splice -uf -k14 -p 0.8 -N 100 "$CWD/genome_index.mmi"  "$SAMPLE_DIR/${sample_id}.fastq" > "$SAMPLE_DIR/${sample_id}.sam"
```

## Quantification

For data derived from long-reads technology the code line is practically the same. Only the long-reads argument has to be activated:
```R
quant <- featureCounts(files = print(paste0(sample_dir,"/",sample_id,".bam")), annot.ext = repeat_gtf,isGTFAnnotationFile = T,
         isLongRead=T, isPairedEnd = FALSE, GTF.attrType = "gene_id",countMultiMappingReads = TRUE,primaryOnly = FALSE, fraction=TRUE,
         strandSpecific = strandness_fc)

```

# Downstream analysis
Differential expression analysis, coexpression analysis, transcriptome assembly and prediction models design steps are common between the three PERREO modes.

## Differential Expression Analysis
Statistical analysis is performed by DESeq2 or edgeR functions using default thresholds: abs(log2FC) > 1 and FDR < 0.05. DESeq2 is a very robust mnethod and is more conservative, while edgeR is specially sensitive for low-expressed genes due to its flexible dispersion estimation. If batch effect = yes, RUVg will be run in case the sample sheet does not contain a batch column indicating the specific cause of this variability in the data. On the other hand, if this column is indicated in the sample sheet, that variable will be directly included in the formula. Different plots and files are exported from this step:<br>
1) Repeat RNAs counts represented with violin plots.<br>
2) PCA (if batch effect is corrected, PCA before and after batch effect reduction are exported).<br>
3) Fold-Change barplots of differentially expressed features (DEFs) for each contrast.<br>
4) Classification of all repeat RNAs identified by repeat classes.<br>
5) Classification of differentially expressed repeat RNAs by repeat classes.<br>
6) Expression matrix is CSV format.<br>
7) Heatmap of DEFs in each contrast.<br>
8) Volcano plot for each contrast.<br>


## Coexpression analysis
WGCNA R package generates coexpression networks from the given expression matrix. Power is automatically assigned when R value is 0.9 and the script is designed to select the three modules with the highest correlation value with respect any of the experimental conditions. The script exports:
1) Information about nodes and edges of the three modules previously mentioned to generate Cytoscape networks.
2) Correlation heatmap between modules and experimental conditions.
3) Correlation heatmap between modules and samples.
4) Repeat RNAs and the module they belong to in txt format.
5) Modules dendrogram.<br>


## Transcriptome assembly

Stringtie2 generates a GTF file for transcriptome assembly of each sample in the study using a combined GTF of genomic and repeat annotations as reference. Strandedness has to be taken into account to perform the assembly. <br>

```bash
stringtie "$BAM" \
        -L \
        -p "$threads" \
        -G "$combined_annotation" \
        -o "$OUTPUT_GTF" \
        --fr \ #rf is strandness=reverse
```

Then, the pipeline performs a preliminar anylsis where the annotations file of each sample is mapped to a GTF containing exonic regions and another GTF containing repetitive elements in order to obtain hybrid transcripts composed of fragments derived from exons and repeats.<br>


## PDF report
The generated report contains the plots previously mentioned in the Differential Expression section.
<img width="6000" height="4500" alt="VolcanoPlot" src="https://github.com/user-attachments/assets/fbcd3692-6819-4f7e-a95f-46432376ff58" />
<img width="2400" height="1800" alt="repetitive_counts_violin_box" src="https://github.com/user-attachments/assets/c7bde717-9a58-43da-89f2-2b35d2c12b94" />
<img width="2400" height="1800" alt="pca" src="https://github.com/user-attachments/assets/a16a83d9-48bf-4048-91bc-e6553f1f734f" />
<img width="2600" height="2000" alt<img width="6000" height="4500" alt="Classification_DEGs" src="https://github.com/user-attachments/assets/94d2dc91-070f-4071-b43b-131f427e96e4" />
<img width="6000" height="4500" alt="Classification_All" src="https://github.com/user-attachments/assets/37aa337a-d45d-4778-be2f-90bb57105c25" />
<img width="6000" height="4500" alt="BarPlotUpDown" src="https://github.com/user-attachments/assets/69d56e63-f1c4-4996-b6a4-8ff4072e35c1" />

