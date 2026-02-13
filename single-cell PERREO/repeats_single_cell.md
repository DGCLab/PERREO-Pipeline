# PERREO version for single-cell data analysis
Single-cell approaches are focused on identifying mRNA. In this way, non-polyadenilated repeats and lncRNAs detection still constitutes an important challenge for this type of experiment. Here we propose a general pipeline to identify repeats RNAs at single-cell resolution.

# Summary
To run this pipeline, it is necessary to provide a preprocessed fastq file per sample. UMIs and barcodes should be previously extracted using specific software as UMI-tools extract. The first step is carried out by STAR aligner, which align the reads with respect to the reference genome applying the two-pass mode. Then, featureCounts is used to tag the corresponding reads and UMI-tools count to perform quantification analysis at UMI-level.

# Downstream analysis



