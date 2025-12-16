#!/usr/bin/env bash
set -euo pipefail

GTF_DIR="$1"
CWD="$2"
genome_gtf="$3"
repeat_gtf="$4"

for gtf in "$GTF_DIR"/*.gtf; do
  sample="$(basename "$gtf" .gtf)"

  # Exones del ensamblado
  awk 'BEGIN{FS=OFS="\t"} $3=="exon"{print}' "$gtf" > tmp.exons.gtf

  # IDs de transcritos con exones que solapan EXONES génicos
  bedtools intersect -u -a tmp.exons.gtf -b <(awk 'BEGIN{FS=OFS="\t"} $3=="exon"{print}' "$CWD/$genome_gtf") \
    | sed -n 's/.*transcript_id "\([^"]*\)".*/\1/p' \
    | sort -u > tmp.tx_gene.ids || true

  # IDs de transcritos con exones que solapan repeticiones
  # (si tu repeat_gtf es todo exon, esto equivale a no filtrar)
  bedtools intersect -u -a tmp.exons.gtf -b <(awk 'BEGIN{FS=OFS="\t"} !/^#/{print}' "$CWD/$repeat_gtf") \
    | sed -n 's/.*transcript_id "\([^"]*\)".*/\1/p' \
    | sort -u > tmp.tx_repeat.ids || true

  # Todos los transcript_id
  awk 'BEGIN{FS=OFS="\t"} $3=="transcript"{print}' "$gtf" \
    | sed -n 's/.*transcript_id "\([^"]*\)".*/\1/p' \
    | sort -u > tmp.tx_all.ids

  total=$(wc -l < tmp.tx_all.ids)

  hybrid=$(comm -12 tmp.tx_gene.ids tmp.tx_repeat.ids | wc -l)
  gene_only=$(comm -23 tmp.tx_gene.ids tmp.tx_repeat.ids | wc -l)
  repeat_only=$(comm -13 tmp.tx_gene.ids tmp.tx_repeat.ids | wc -l)

  frac=$(awk -v h="$hybrid" -v t="$total" 'BEGIN{if(t>0) printf "%.4f", h/t; else print 0}')

done