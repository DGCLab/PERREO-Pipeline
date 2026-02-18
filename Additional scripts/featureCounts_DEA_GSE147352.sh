/usr/bin/time -v ./featurecounts_dea.sh -sample_list samplesheet.txt -reference_genome nanopore_genome.fa -genome_gtf genome_reference_nanopore_annotation_2.gtf \
  -repeat_gtf repeatmasker_nanopore_annotation_3.gtf -threads 14 -trimming simple -adapt_r1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -adapt_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -project_name "repeat_rna_glioblastoma_project" -remove_duplicates false -batch no -method DESeq2 -FDR 0.01 -positive_class "GBM" > PERREO_performance.log 2>&1
