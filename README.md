# PERREO-Pipeline

Instructions:<br> 

Create a conda environment

```bash
conda create --name repeat_rnaseq
conda activate repeat_rnaseq

```

Software requirements to run this pipeline for bash and R are the following: <br> 

```bash
conda install bioconda::multiqc
conda install bioconda::fastqc
conda install bioconda::star
conda install bioconda::picard
conda install bioconda::cutadapt
conda install bioconda::gatk4
conda install conda-forge::natsort
conda install conda-forge::r-base
conda install conda-forge::parallel
conda install conda-forge::r-rjava
conda install bioconda::nanoplot
conda install bioconda::minimap2
conda install bioconda::nanocount


```

```bash
Rscript -e 'BiocManager::install("Rsubread")'
Rscript -e 'BiocManager::install("DESeq2")'
Rscript -e 'BiocManager::install("edgeR")'
Rscript -e 'BiocManager::install("rtracklayer")'
Rscript -e 'BiocManager::install("EDASeq")'
Rscript -e 'BiocManager::install("RUVSeq")'
Rscript -e 'BiocManager::install("WGCNA")'
Rscript -e 'install.packages("tidyverse",repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("readxl",repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("openxlsx", repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("dplyr", repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("readr", repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("tidyverse", repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("ggplot2", repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("pheatmap", repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("ggpubr", repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("ggrepel", repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("mlbench", repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("caret", repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("caretEnsemble", repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("randomForest", repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("glmnet", repos="https://cloud.r-project.org")'

```

Repeat annotations is a complex task and it can take many hours. In this case, we provide a worflow to obtain repeat annotations using RepeatModeler, MCHelper and Repeatmasker.<br> 
<br>
In case the user prefers to use prebuilt annotations, our pipeline can be run providing customized annotations publicly available in different platforms.<br> 

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

To run this pipeline it is necessary to indicate different parameters depending on the selected mode:<br> 

Paired-end short-reads RNA-seq mode<br> 

```bash
-sample_list
-reference_genome
-genome_gtf
-repeat_gtf: must contain a "repeat_class" column in order to study the type of repetitive elements identified in the analysis.
-threads #(default: threads=8)
-adapt_r1: If reads are already trimmed, the user can ignore this argument and trimming will be also performed to obtain high quality reads.
-adapt_r2: If reads are already trimmed, the user can ignore this argument and trimming will be also performed to obtain high quality reads.
-trimming_type: "trimming_extra" if  "trimming_simple". #trimming_simple/trimming_extra (default: trimming_type=trimming_simple)
-project_name
-remove_duplicates #yes/no
-batch_effect #yes/no
-method #edgeR/DESeq2
-k_num #(default: k_num=2)
```

Paired-end short-reads RNA-seq mode<br> 
```bash

```

Long-reads direct RNA-seq mode<br> 
```bash

```

The samplesheet structure should be like this:<br>
```bash
sample	strandedness	condition
SRR14506659	reverse	ESO
SRR14506660	reverse	ESO
SRR14506661	reverse	ESO
SRR14506662	reverse	ESO
SRR14506859	reverse	HC
SRR14506860	reverse	HC
SRR14506861	reverse	HC
SRR14506862	reverse	HC
```
*** In the case the data belong to more than two experimental conditions and batch effect cause is known, another column called "batch" must be included in the samplesheet indicating the origin batch. <br>
<br>
For example: If samples come from different hospital, this column should be included and the hospital of origin has to be indicated for each sample.<br>
<br>
The pipeline can be run on both single-end and paired-end data. There are to main scripts to run the analysis depending on the type of data the user has to analyze.<br>
<br>
-Paired-end data pipeline: repRPI_PE.sh<br>
-Single-end data pipeline: repRPI_SE.sh<br>
<br>
Regarding the reference genome, the user can decide which reference genome and which annotations to use as the software does not provide any of these files. Consequently, this pipeline is applicable to any organism whose genome is sequenced and annotated.<br> 
For human experiments, although reference and annotations can be freely chosen, we recommend to use T2T genome and its corresponding annotation as it is the most complete human reference, where repetitive regions are better described.<br>

Here we provide the link to the T2T fasta file and GTF annotations to use in this workflow:https://github.com/marbl/CHM13.<br>

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

![FIGURA PAPER 1](https://github.com/user-attachments/assets/eedee39c-b058-4040-b282-a394f08500ec)


# Trimming
In this step there are two main options: simple trimming with cutadapt and a more complex trimming performed firstly with cutadapt and then with trimGC.py script in order to remove additional GC nucleotides added by specific sequencing kits. <br>
<br>
Another possibility is that adapters have been removed previously and the trimming step should only be taken into account to remove low quality reads. In this case, user must not include -adapt_r1 and -adapt_r2 arguments while running the pipeline.<br>


# Alignment
STAR aligner parameters are indicated by default allowing the multimapping and removing reads with more than a 10% of mismatches. outFilterMismatchNoverLmax parameter can be changed depending on the user experimental design and conditions.<br>

```bash
STAR --runThreadN $threads \
      --genomeDir "$GENOME_DIR" \
      --readFilesIn "$trimmed1" "$trimmed2" \
      --outFileNamePrefix "$MAP_DIR/${sample_id}_" \
      --outSAMtype BAM SortedByCoordinate \
      --outFilterMultimapNmax 500 \
      --winAnchorMultimapNmax 500 \
      --outFilterMismatchNoverLmax 0.1 \
      --outSAMmultNmax 500 \
      --outMultimapperOrder Random \
      --runRNGseed 42 \
      --outSAMattributes NH HI AS nM NM MD
```

# Duplications analysis
Duplications removal is not recommended generally in RNA-seq data analysis. However, there are situations where the percentage of duplications is too high and remove them is a real option.


# Quantification


# Differential Expression Analysis
<img width="6000" height="4500" alt="VolcanoPlot" src="https://github.com/user-attachments/assets/fbcd3692-6819-4f7e-a95f-46432376ff58" />
<img width="2400" height="1800" alt="repetitive_counts_violin_box" src="https://github.com/user-attachments/assets/c7bde717-9a58-43da-89f2-2b35d2c12b94" />
<img width="2400" height="1800" alt="pca" src="https://github.com/user-attachments/assets/a16a83d9-48bf-4048-91bc-e6553f1f734f" />
<img width="2600" height="2000" alt<img width="6000" height="4500" alt="Classification_DEGs" src="https://github.com/user-attachments/assets/94d2dc91-070f-4071-b43b-131f427e96e4" />
<img width="6000" height="4500" alt="Classification_All" src="https://github.com/user-attachments/assets/37aa337a-d45d-4778-be2f-90bb57105c25" />
<img width="6000" height="4500" alt="BarPlotUpDown" src="https://github.com/user-attachments/assets/69d56e63-f1c4-4996-b6a4-8ff4072e35c1" />

