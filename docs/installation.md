# Installation

PERREO requires a conda environment with the command-line tools and R packages used across the pipeline. Create the environment from the provided YAML file:

```bash
conda env create -f perreo.yml
```

Activate the environment before running the pipeline:

```bash
conda activate perreo
```

The repository includes the main interface, mode-specific entry scripts, and supporting R and shell scripts:

```text
app.R
perreo.yml
short_reads_PE/
short_reads_SE/
long_reads/
Additional scripts/
```

## External Inputs

PERREO does not provide reference genomes or annotation files. Users must provide:

- A reference genome in FASTA format.
- A genome annotation GTF for short-read modes.
- A repeat annotation GTF with a `repeat_class` attribute.

For human experiments, the T2T reference and annotations from the CHM13 project are recommended because repetitive regions are better represented.

## Preanalysis

Before launching PERREO, run FastQC and MultiQC on raw FASTQ files. These reports help decide trimming settings and whether duplicate removal is appropriate.

```bash
mkdir -p QC/fastqc_raw

find . -type f \( -iname "*.fastq" -o -iname "*.fq" -o -iname "*.fastq.gz" -o -iname "*.fq.gz" \) -print0 \
| xargs -0 -n 1 -P 4 fastqc -t 2 -o QC/fastqc_raw

multiqc QC/fastqc_raw
```
