# Downstream Analysis

The downstream workflow is shared across SR-SE, SR-PE, and LR modes, with mode-specific quantification settings.

## Differential Expression

PERREO runs differential expression analysis with either DESeq2 or edgeR. Default thresholds are:

- `abs(log2FC) > 1`
- `FDR < 0.05`

If batch correction is enabled and the samplesheet contains a `batch` column, the batch term is included in the model design. If batch correction is enabled without a `batch` column, PERREO applies RUVg to reduce unwanted variation.

The differential expression step exports:

1. Repeat RNA count violin plots.
2. PCA plots, including before and after correction when batch correction is applied.
3. Fold-change bar plots for differentially expressed features.
4. Repeat RNA classification by repeat class.
5. Differentially expressed repeat RNA classification by repeat class.
6. Expression matrices in CSV or TXT format.
7. Heatmaps of differentially expressed features.
8. Volcano plots for each contrast.

## Coexpression

WGCNA is used to build coexpression networks from the expression matrix. The scripts automatically select the soft-thresholding power and export the three modules with the strongest correlation to experimental conditions.

Outputs include node and edge tables for Cytoscape, module-condition heatmaps, module-sample heatmaps, repeat RNA module assignments, and module dendrograms.

## Transcriptome Assembly

StringTie2 assembles transcripts using combined genome and repeat annotations. Strandedness is passed to the assembly step.

```bash
stringtie "$BAM" \
  -L \
  -p "$threads" \
  -G "$combined_annotation" \
  -o "$OUTPUT_GTF"
```

The pipeline then identifies candidate hybrid transcripts by intersecting sample-level transcript annotations with exonic and repetitive regions.

## Prediction Models

PERREO can train GLMnet and Random Forest models to evaluate whether repeat RNA profiles classify experimental conditions. Prediction models require enough samples per class and a defined positive class.

A common error in very small or imbalanced datasets is:

```text
Error in roc.default(...): No case observation.
```

If this persists, review the number of samples per class, consider merging underrepresented classes, or exclude selected samples only from the prediction model step after upstream analyses are complete.
