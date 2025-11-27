#!/bin/bash
# Installation of required software
# Conda installation
conda install bioconda::multiqc
conda install bioconda::fastqc
conda install bioconda::star
conda install bioconda::picard
conda install bioconda::cutadapt
conda install bioconda::gatk4
conda install conda-forge::natsort
conda install conda-forge::r-base
conda install conda-forge::parallel
conda install -c conda-forge r-rjava
conda install bioconda::nanoplot
conda install bioconda::minimap2
conda install bioconda::nanocount

# Installation with R
#Rscript -e 'BiocManager::install("Rsubread")'
#Rscript -e 'BiocManager::install("DESeq2")'
#Rscript -e 'BiocManager::install("edgeR")'
#Rscript -e 'BiocManager::install("rtracklayer")'
#Rscript -e 'BiocManager::install("EDASeq")'
#Rscript -e 'BiocManager::install("RUVSeq")'
#Rscript -e 'BiocManager::install("WGCNA")'
#Rscript -e 'install.packages("tidyverse",repos="https://cloud.r-project.org")'
#Rscript -e 'install.packages("readxl",repos="https://cloud.r-project.org")'
#Rscript -e 'install.packages("openxlsx", repos="https://cloud.r-project.org")'
#Rscript -e 'install.packages("dplyr", repos="https://cloud.r-project.org")'
#Rscript -e 'install.packages("readr", repos="https://cloud.r-project.org")'
#Rscript -e 'install.packages("tidyverse", repos="https://cloud.r-project.org")'
#Rscript -e 'install.packages("ggplot2", repos="https://cloud.r-project.org")'
#Rscript -e 'install.packages("pheatmap", repos="https://cloud.r-project.org")'
#Rscript -e 'install.packages("ggpubr", repos="https://cloud.r-project.org")'
#Rscript -e 'install.packages("ggrepel", repos="https://cloud.r-project.org")'
#Rscript -e 'install.packages("mlbench", repos="https://cloud.r-project.org")'
#Rscript -e 'install.packages("caret", repos="https://cloud.r-project.org")'
#Rscript -e 'install.packages("caretEnsemble", repos="https://cloud.r-project.org")'
#Rscript -e 'install.packages("randomForest", repos="https://cloud.r-project.org")'
#Rscript -e 'install.packages("glmnet", repos="https://cloud.r-project.org")'

conda install bioconda::bioconductor-rsubread
conda install bioconda::bioconductor-DESeq2
conda install bioconda::bioconductor-edgeR
conda install bioconda::bioconductor-rtracklayer
conda install bioconda::bioconductor-EDASeq
conda install bioconda::bioconductor-RUVSeq
conda install bioconda::r-wgcna
Rscript -e 'install.packages("tidyverse", repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("ggplot2", repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("pheatmap", repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("ggpubr", repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("ggrepel", repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("mlbench", repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("caret", repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("caretEnsemble", repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("randomForest", repos="https://cloud.r-project.org")'
conda install -c conda-forge r-glmnet

