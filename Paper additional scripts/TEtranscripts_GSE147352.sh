/usr/bin/time -v TEtranscripts -t $(xargs -a GBM.txt) -c $(xargs -a HC.txt) \
  --GTF genome_reference_nanopore_annotation_2.gtf --TE T2T_CHM13_v2_rmsk_TE.gtf --stranded reverse --mode multi \
  --project GLIO_TEtranscripts --outdir Results --padj 0.01 --foldchange 1 --sortByPos > tetranscripts_log 2>&1
