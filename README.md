# PERREO-Pipeline

# Software requirements:<br> 

It is necessary to create a conda environment with all the required packages and programs to perform the whole analysis. To do it, the user can run this yaml file. In this way, the complete environment will be installed automatically. 
```bash
conda env create -f perreo.yml

```

# Repeat annotation

Firstly, it is necessary to build a database from a reference genome and to run RepeatModeler using the following code lines:<br> 

```bash
BuildDatabase -name homo_sapiens reference_genome.fa

RepeatModeler -database ../repeat_annotations/homo_sapiens \
                   -threads ${SLURM_CPUS_PER_TASK:-24} \
                   -LTRStruct
```

Then, it is recommended to curate the RepeatModeler outputs, and we propose MCHelper as software to curate the annotations automatically.  To run the software properly it is important to download and prepare all the needed libraries as described in this link: https://github.com/GonzalezLab/MCHelper?tab=readme-ov-file#linuxwindows.<br> 

```bash
python MCHelper.py -l consensi.fa.classified -o MCHelper_out -g nanopore_genome.fa --input_type fasta -b db/mammalia_odb10

```

Finally...........<br> 


# Data preparation<br> 

The suitable folders structure for the correct performance of the software should be the following:<br> 

```bash
PROJECT_FOLDER/
├─ SAMPLES
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
For human experiments, although reference and annotations can be freely chosen, we recommend to use T2T genome and its corresponding annotation as it is the most complete human reference, where repetitive regions are better described.<br>

Here we provide a link to the T2T fasta file and GTF annotations to use in this workflow:https://github.com/marbl/CHM13.<br>

In the case you need to download other fasta files and specific GTF annotations of T2T, you can also use the UCSC table browser platform: https://genome.ucsc.edu/cgi-bin/hgTables.<br>
<br>
# Required preanalysis<br>

Before running the pipeline you should perform fastqc and multiqc to have the enough information to make some decisions with respect to quality metrics before trimming the reads and selecting or not the remove duplicates option of MarkDuplicates.<br>


To run multiqc on fastqc outputs, the user can use these command lines:

```bash
mkdir -p QC/fastqc_raw

# Process .fastq, .fq, .fastq.gz y .fq.gz
find . -type f \( -iname "*.fastq" -o -iname "*.fq" -o -iname "*.fastq.gz" -o -iname "*.fq.gz" \) -print0 \
| xargs -0 -n 1 -P 4 fastqc -t 2 -o QC/fastqc_raw
```

# Workflow summary<br>

The flow diagram describes the different steps taken into account in this pipeline.<br>
![workflow_perreo](https://github.com/user-attachments/assets/a737c6e9-a09d-457a-86ca-304b3e702e7f)


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
| Trimming type     |  trimming_simple (not needed for PERREO LR)     | 
| Mismatch align     | 0.05     | 
| Trimming quality threshold     | 30     | 
| Minimum reads length (for trimming)     | 16     | 
| Maximum reads length (for trimming)     |      | 
| K_num     | 2     | 
| log2FC     | 1     | 
| FDR     | 0.05     | 
| Prediction model     | yes     | 


For SR-LR mode, the parameters trimming_type, trimming_quality_threshold, min_length_trim and max_length_trim are not considered as trimming for long-reads should be performed during the basecalling before obtaining the fastq files.


# PERREO SR-PE and SR-SE
The required arguments for both paired-end-derived data are the following:<br> 

For SR-PE:<br> 
```text
-sample_list                  Sample sheet with sample, strandedness, condition and batch (if necessary).
-reference_genome             Genome file in fasta.
-genome_gtf                   Genome annotations in GTF format.
-repeat_gtf                   Repeat annotations in GTF format. It must contain a "repeat_class" column in order to study the type
                              of repetitive elements identified in the analysis.
-threads                      Number of threads used for the process (default: 8)
-adapt_r1                     Adapter 1 sequence. If reads are already trimmed, the user can ignore this argument and trimming will
                              be also performed to obtain high quality reads.
-adapt_r2                     Adapter 2 sequence. If reads are already trimmed, the user can ignore this argument and trimming will
                              be also performed to obtain high quality reads.
-trimming_type                trimming_simple/trimming_extra. The second must be selected if it is know that the kit used adds extra
                              GC nucleotides (default: trimming_simple).
-trimming_quality_threshold   Minimum quality permitted for reads to be kept after trimming (default: 30).
-min_length_trim              Minimum reads length to not discard them after trimming (default: 16).
-max_length_trim              Maximum reads length permitted for reads to be kept after trimming (default: ).
-mismatch_align               Percentage of mismatches permitted for reads to be kept during alignment (default: 0.05).
-project_name
-remove_duplicates            yes/no. Generally it is not recommended to discard duplicates unless the proportion is really high.
-batch_effect                 yes/no. If yes, and batch column is included in the sample sheet, the batch effect will be reduced
                              based on that column. If yes but batch column does not exist, RUVg will reduce the undesired variability.
-method                       edgeR/DESeq2. 
-k_num                        K number for RUVg-based batch effect reduction. As far as the dataset size increases, it is recommended
                              to also increase the k number (default: 2).
-log2FC                       Log2-transformed fold-change threshold for Differential Expression Analysis (default: 1).
-FDR                          Adjusted p-value threshold for Differential Expression Analysis (default: 0.05).
-prediction_model             yes/no (default: yes). When activated, it only will design prediction models in case the number of samples
                              is higher than 40. 
```
For SR-SE:<br> 

```bash

```


## Trimming
In this step there are two main options: simple trimming with cutadapt and a more complex trimming performed firstly with cutadapt and then with trimGC.py script in order to remove additional GC nucleotides added by specific sequencing kits. <br>
<br>
If adapters had been removed previously and the trimming step should only be taken into account to remove low quality reads, user must not include -adapt_r1 and -adapt_r2 arguments while running the pipeline.<br>


## Alignment
STAR aligner parameters are indicated by default allowing the multimapping and removing reads with more than a 5% of mismatches. This parameter can be changed depending on the user experimental design and conditions.<br>

```bash
STAR --runThreadN $threads \
      --genomeDir "$GENOME_DIR" \
      --readFilesIn "$trimmed1" "$trimmed2" \
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

## Duplications analysis
Duplications removal is not recommended generally in RNA-seq data analysis. However, there are situations where the percentage of duplications is too high and remove them is a real option.


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
-prediction_model      yes/no (default: yes). When activated, it only will design prediction models in case the number of samples
                       is higher than 40. 
```

In this case, trimming parameters are not required as adapters and barcodes removal is carried out during the basecalling step and quality control is carried out using NanoPlot.<br>

## Alignment
Long reads are aligned against the reference genome using minimap2 with -ax splice -uf and -k14 parameters. Multimapping is also allowed indicating -N 100 parameter. 

```bash
minimap2 -t 14 -ax splice -uf -k14 -p 0.8 -N 100 "$CWD/genome_index.mmi"  "$SAMPLE_DIR/${sample_id}.fastq" > "$SAMPLE_DIR/${sample_id}.sam"
```

# Downstream analysis
Quantification, differential expression analysis, coexpression analysis, transcriptome assembly and prediction models design steps are common between the three PERREO modes.

## Quantification
FeatureCounts performs features quantification allowing multimapping reads counts and fraction. It uses the strandedness indicated in the sample sheet.<br>

The following code line is run for data obtained from short reads sequencing techonology:

```R
quant <- featureCounts(files = print(paste0(MAP_DIR,"/",sample_id,"_marked_duplicates_STAR.bam")), annot.ext = repeat_gtf,
         isGTFAnnotationFile = T, isPairedEnd = TRUE, GTF.attrType = "gene_id",countMultiMappingReads = TRUE,primaryOnly = FALSE,
         fraction=TRUE,strandSpecific = strandness_fc)
```

For data derived from long-reads technology the code line is practically the same. Only the long-reads argument has to be activated:
```R
quant <- featureCounts(files = print(paste0(sample_dir,"/",sample_id,".bam")), annot.ext = repeat_gtf,isGTFAnnotationFile = T,
         isLongRead=T, isPairedEnd = FALSE, GTF.attrType = "gene_id",countMultiMappingReads = TRUE,primaryOnly = FALSE, fraction=TRUE,
         strandSpecific = strandness_fc)

```

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

## Transcriptome assembly
Stringtie2 generates a GTF file for transcriptome assembly of each sample in the study. Strandedness has to be taken into account to perform the assembly. <br>

```bash
stringtie "$BAM" \
        -L \
        -p "$threads" \
        -G "$repeat_annotation" \
        -o "$OUTPUT_GTF" \
        --fr \ #rf is strandness=reverse
```

## Coexpression analysis
WGCNA R package generates coexpression networks from the given expression matrix. Power is automatically assigned when R value is 0.9 and the script is designed to select the three modules with the highest correlation value with respect any of the experimental conditions. The script exports:
1) Information about nodes and edges of the three modules previously mentioned to generate Cytoscape networks.
2) Correlation heatmap between modules and experimental conditions.
3) Correlation heatmap between modules and samples.
4) Repeat RNAs and the module they belong to in txt format.
5) Modules dendrogram.  


## PDF report
The generated report contains the plots previously mentioned in the Differential Expression and Coexpression analysis sections.
<img width="6000" height="4500" alt="VolcanoPlot" src="https://github.com/user-attachments/assets/fbcd3692-6819-4f7e-a95f-46432376ff58" />
<img width="2400" height="1800" alt="repetitive_counts_violin_box" src="https://github.com/user-attachments/assets/c7bde717-9a58-43da-89f2-2b35d2c12b94" />
<img width="2400" height="1800" alt="pca" src="https://github.com/user-attachments/assets/a16a83d9-48bf-4048-91bc-e6553f1f734f" />
<img width="2600" height="2000" alt<img width="6000" height="4500" alt="Classification_DEGs" src="https://github.com/user-attachments/assets/94d2dc91-070f-4071-b43b-131f427e96e4" />
<img width="6000" height="4500" alt="Classification_All" src="https://github.com/user-attachments/assets/37aa337a-d45d-4778-be2f-90bb57105c25" />
<img width="6000" height="4500" alt="BarPlotUpDown" src="https://github.com/user-attachments/assets/69d56e63-f1c4-4996-b6a4-8ff4072e35c1" />

