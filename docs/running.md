# Running PERREO

PERREO can be launched through the Shiny graphical interface or through the mode-specific shell scripts.

## Graphical Interface

Place `app.R` in the project directory and run:

```bash
export PATH="/home/user/project_directory:$PATH"
Rscript app.R
```

The interface opens in a web browser and allows users to select the sequencing mode, upload files, set parameters, and monitor execution messages in the terminal.

## Command-Line Entry Points

The repository provides one entry script per mode:

| Mode | Script |
| --- | --- |
| SR-PE | `short_reads_PE/perreo_srpe.sh` |
| SR-SE | `short_reads_SE/perreo_srse.sh` |
| LR | `long_reads/perreo_lr.sh` |

Each entry script expects to run from the project directory containing `samples/`, the samplesheet, reference files, and the mode-specific `scripts_*` directory.

## Common Arguments

| Argument | Description | Default |
| --- | --- | --- |
| `-sample_list` | Tab-separated samplesheet. | Required |
| `-reference_genome` | Reference genome FASTA. | Required |
| `-repeat_gtf` | Repeat annotation GTF with `repeat_class`. | Required |
| `-threads` | Number of processing threads. | `8` |
| `-project_name` | Project name used in outputs and reports. | Required |
| `-method` | Differential expression method: `edgeR` or `DESeq2`. | Required |
| `-log2FC` | Absolute log2 fold-change threshold. | `1` |
| `-FDR` | Adjusted p-value threshold. | `0.05` |
| `-prediction_model` | Whether to run prediction models: `yes` or `no`. | `no` |
| `-positive_class` | Positive class for prediction models. | Required if prediction models are enabled |

Short-read modes use `-batch`; long-read mode uses `-batch_effect`.
