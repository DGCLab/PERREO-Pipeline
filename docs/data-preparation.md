# Data Preparation

PERREO expects a project directory with sample folders, reference files, annotations, the samplesheet, and the scripts for the selected mode.

```text
PROJECT_FOLDER/
├── samples/
│   ├── sample1/
│   └── sample2/
├── reference_genome.fa
├── genome.gtf
├── repeats.gtf
├── samplesheet.txt
└── scripts/
```

The mode-specific script folder from this repository must be available from the directory where the pipeline is executed.

## Samplesheet

The samplesheet is a tab-separated file with at least three columns:

| sample | strandedness | condition |
| --- | --- | --- |
| SRR14506659 | reverse | ESO |
| SRR14506660 | reverse | ESO |
| SRR14506859 | reverse | HC |

Use `forward`, `reverse`, or `unstranded` in the `strandedness` column. If strandedness is unknown, infer it before using extra trimming, because `trimGC.py` depends on the library orientation.

## Batch Column

If samples come from more than two experimental conditions and a known batch source exists, add a `batch` column:

| sample | strandedness | condition | batch |
| --- | --- | --- | --- |
| SRR14506659 | reverse | ESO | Hosp_A |
| SRR14506660 | reverse | ESO | Hosp_B |
| SRR14506859 | reverse | HC | Hosp_C |

When batch correction is enabled, PERREO uses this column if present. If batch correction is enabled without a `batch` column, RUVg is used to reduce unwanted variation.

## Annotation Consistency

Genome and repeat annotations must use matching chromosome names. For example, avoid mixing `chrX` in one file with `X` in another file. Standardize chromosome names before running the pipeline.
