# Outputs

PERREO writes intermediate and final files into sample-level folders and a project-level `Results/` directory.

## Sample-Level Outputs

Each sample may include:

- Trimmed FASTQ files.
- Alignment files in BAM format.
- Alignment indexes.
- Quantification tables.
- Quality control reports.
- Transcriptome assembly files.

## Project-Level Outputs

The `Results/` directory contains combined outputs from downstream analyses, including:

- Count matrices.
- Differential expression result tables.
- PCA plots.
- Volcano plots.
- Fold-change bar plots.
- Repeat class summary plots.
- Heatmaps.
- WGCNA module files.
- Cytoscape node and edge tables.
- Hybrid transcript summaries.
- Prediction model metrics and plots when enabled.
- PDF reports.

## PDF Report

The report includes the main plots from differential expression analysis:

<img width="6000" height="4500" alt="Volcano plot" src="https://github.com/user-attachments/assets/fbcd3692-6819-4f7e-a95f-46432376ff58" />

<img width="2400" height="1800" alt="Repeat count violin plot" src="https://github.com/user-attachments/assets/c7bde717-9a58-43da-89f2-2b35d2c12b94" />

<img width="2400" height="1800" alt="PCA plot" src="https://github.com/user-attachments/assets/a16a83d9-48bf-4048-91bc-e6553f1f734f" />

<img width="6000" height="4500" alt="Repeat classification plot" src="https://github.com/user-attachments/assets/37aa337a-d45d-4778-be2f-90bb57105c25" />

<img width="6000" height="4500" alt="Up/down bar plot" src="https://github.com/user-attachments/assets/69d56e63-f1c4-4996-b6a4-8ff4072e35c1" />
