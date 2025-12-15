#!/usr/bin/env bash
#set -euo pipefail

threads="$1"
PREPDE_SCRIPT="$2"
CWD="$3"
genome_gtf="$4"
repeat_gtf="$5"
sample_list="$6"
sample_list="$CWD/$sample_list"
strandedness="$7"


# Folder where all the GTFs from stringtie are stored
GTF_DIR="$CWD/Transcriptome_assembly"                       

# Genome GTF
REF_GTF="$CWD/$genome_gtf"               

# SAMPLES directory
SAMPLES_DIR="$CWD/SAMPLES"

# Repeats annotations
RM_GTF="$CWD/$repat_gtf"

# Outputs
SECOND_DIR="$CWD/Transcriptome_assembly_novels"
MERGELIST="$SECOND_DIR/mergelist.txt"
MERGED_GTF="$SECOND_DIR/merged.gtf"

# prepDE.py script location
PREPDE_SCRIPT="$PREPDE_SCRIPT"

#Creating merge list from the GTFs stored in Transcriptome_assembly folder

ls "$GTF_DIR"/*.gtf > "$MERGELIST"


# 1) Performing stringtie --merge on all the GTFs

echo "[INFO] Ejecutando stringtie --merge..."
stringtie --merge \
  -G "$REF_GTF" \
  -o "$MERGED_GTF" \
  "$MERGELIST"

echo "[OK] merged.gtf generado: $MERGED_GTF"


# 2) gffcompare for CLASS_CODE (NOVEL vs GENES)
# This step is necessary to mark novels transcripts with respect to genome annotations

gffcompare -r "$REF_GTF" -o gffcmp "$MERGED_GTF"


# 3) New stringtie analysis using the merge gtf

mkdir -p "$SECOND_DIR"
echo "sample,sample_gff" > sample_list.csv

# sample_list: to take the first column with sample_ids
tail -n +2 "$sample_list" | while IFS=$'\t' read -r sample condition bam_col strandedness; do
    sample_dir="$SAMPLES_DIR/$sample"
    [ -d "$sample_dir" ] || continue

    echo "[INFO] Procesando muestra $sample"

    bam="$sample_dir/alignment/${sample}.bam"
    if [ ! -f "$bam" ]; then
        echo "[WARN] Expected BAM not found: $bam. Sample skipped."
        continue
    fi
    echo "    BAM: $bam"

    out_gtf="$SECOND_DIR/${sample}.gtf"

      if [ "$strandedness" = "forward" ]; then
    echo "Using forward parameter"
        stringtie "$bam" \
      -G "$MERGED_GTF" \
      -e -B -p "$threads" \
      -o "$out_gtf" --fr
    fi

    if [ "$strandedness" = "reverse" ]; then
    echo "Using reverse parameter"
        stringtie "$bam" \
      -G "$MERGED_GTF" \
      -e -B -p "$threads" \
      -o "$out_gtf" --rf
    fi



    # Add to sample_list.csv for prepDE.py
    echo "${sample},${out_gtf}" >> sample_list.csv
done


# 4) Generating count matrixes with prepDE.py3

echo "[INFO] Ejecutando prepDE.py..."
python "$PREPDE_SCRIPT" -i sample_list.csv

echo "[OK] gene_count_matrix.csv and transcript_count_matrix.csv generated"


# 5) Overlapping transcripts-TE using only GTF

# 5.1. Exons of merged.gtf (observed transcripts)
echo "[INFO] Extracting exons from merged.gtf..."
awk '$3=="exon"' "$MERGED_GTF" > merged_exons.gtf

# 5.2. Exons of repeatmasker.gtf (repeats features)
echo "[INFO] Extracting exons from repeats.gtf..."
awk '$3=="exon"' "$RM_GTF" > repeat_exons.gtf

# 5.3. bedtools intersect directamente sobre GTF (-wa -wb para conservar atributos)
echo "[INFO] Intersecting transcripts and repeats exons..."
bedtools intersect -wa -wb \
  -a merged_exons.gtf \
  -b repeat_exons.gtf \
  > TE_exon_overlap.gtf


echo "[OK] TE_exon_overlap.gtf generated (GTF + GTF with overlapping information"

cat << 'EOF'
[SUMMARY]
- merged.gtf                 -> catálogo global de transcritos
- gffcmp.annotated.gtf       -> class_code (novel vs genes)
- ${SECOND_DIR}/*.gtf        -> 2ª pasada de StringTie (un GTF por muestra)
- gene_count_matrix.csv
- transcript_count_matrix.csv
- merged_exons.gtf           -> exones de transcritos (StringTie)
- repeat_exons.gtf           -> exones de TE
- TE_exon_overlap.gtf        -> lineas exon_tx + exon_TE solapados
EOF
