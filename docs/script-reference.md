# Script Reference

This page summarizes the repository layout and the main scripts used by each PERREO mode.

## Top-Level Files

| File | Purpose |
| --- | --- |
| `README.md` | Repository overview and original usage notes. |
| `app.R` | Shiny graphical interface for launching PERREO. |
| `perreo.yml` | Conda environment definition. |
| `consensus_search.py` | Utility for consensus sequence window analysis. |
| `logging.sh` | Shared logging helper. |

## Mode Entry Scripts

| Mode | Entry script |
| --- | --- |
| SR-PE | `short_reads_PE/perreo_srpe.sh` |
| SR-SE | `short_reads_SE/perreo_srse.sh` |
| LR | `long_reads/perreo_lr.sh` |

## Short-Read Supporting Scripts

The paired-end and single-end modes contain similar supporting scripts in `scripts_PE/` and `scripts_SE/`:

| Script | Purpose |
| --- | --- |
| `script_trimming_*.sh` | cutadapt-based trimming. |
| `script_map_*.sh` | STAR alignment. |
| `script_markduplicates.sh` | Optional duplicate marking/removal. |
| `quant.R` | featureCounts quantification. |
| `merge_quant.R` | Count matrix merge. |
| `dea.R` | Two-condition differential expression analysis. |
| `dea_multicond.R` | Multi-condition differential expression analysis. |
| `WGCNA.R` | Coexpression analysis. |
| `stringtie2.sh` | Transcriptome assembly. |
| `prediction_model.R` | GLMnet and Random Forest prediction models. |

## Long-Read Supporting Scripts

Long-read scripts live in `long_reads/scripts_long_reads/` and mirror the downstream short-read workflow while using minimap2 and long-read featureCounts settings for upstream processing.

## Additional Scripts

The `Additional scripts/` directory contains analysis examples and auxiliary workflows for external datasets, including TEtranscripts, Salmon, featureCounts, and differential expression examples.
