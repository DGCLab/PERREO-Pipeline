#!/usr/bin/env bash
#set -euo pipefail

threads="$1"
PREPDE_SCRIPT="$2"
CWD="$3"
genome_gtf="$4"
repeat_gtf="$5"
sample_list="$6"
sample_list="$CWD/$sample_list"


##############################
# PARÁMETROS A EDITAR
##############################

# Carpeta donde tienes TODOS los GTF de la 1ª pasada de StringTie
GTF_DIR="$CWD/Transcriptome_assembly"                       # <-- tu carpeta con los .gtf de cada muestra

# Anotación de genes SIN TEs (Gencode/Ensembl)
REF_GTF="$CWD/$genome_gtf"               # pon aquí tu GTF de genes

# Directorio con las muestras y sus BAM
# Estructura esperada:
# SAMPLES/<sample>/alignment/<sample>_Aligned.sortedByCoord.out.bam
SAMPLES_DIR="$CWD/SAMPLES"

# Fichero RepeatMasker en BED (mismo genoma)
RM_GTF="$CWD/$repat_gtf"

# Salidas
MERGELIST="mergelist.txt"
MERGED_GTF="merged.gtf"
SECOND_DIR="$CWD/Transcriptome_assembly_novels"

# Hilos
THREADS="$threads"

# Ruta al prepDE.py
PREPDE_SCRIPT="$PREPDE_SCRIPT"

##############################
# 0) CREAR MERGELIST.TXT DESDE GTF_DIR
##############################

echo "[INFO] Creando mergelist.txt desde $GTF_DIR..."
ls "$GTF_DIR"/*.gtf > "$MERGELIST"
echo "[OK] mergelist.txt con $(wc -l < "$MERGELIST") GTF"

##############################
# 1) MERGE DE TRANSCRITOS
##############################

echo "[INFO] Ejecutando stringtie --merge..."
stringtie --merge \
  -G "$REF_GTF" \
  -o "$MERGED_GTF" \
  "$MERGELIST"

echo "[OK] merged.gtf generado: $MERGED_GTF"

##############################
# 2) gffcompare PARA CLASS_CODE (NOVEL vs GENES)
##############################

echo "[INFO] Ejecutando gffcompare..."
gffcompare -r "$REF_GTF" -o gffcmp "$MERGED_GTF"

echo "[OK] gffcmp.annotated.gtf generado (class_code, ref_gene_id, etc.)"

##############################
# 3) 2ª PASADA DE STRINGTIE CON merged.gtf
#    (TODOS LOS GTF DIRECTAMENTE EN SECOND_DIR)
##############################

mkdir -p "$SECOND_DIR"
echo "sample,sample_gff" > sample_list.csv

# sample_list: fichero con cabecera en la 1ª línea (p.ej. sample<TAB>condition<TAB>bam)
# aquí solo usamos la 1ª columna (sample)
tail -n +2 "$sample_list" | while IFS=$'\t' read -r sample condition bam_col; do
    sample_dir="$SAMPLES_DIR/$sample"
    [ -d "$sample_dir" ] || continue

    echo "[INFO] Procesando muestra $sample"

    bam="$sample_dir/alignment/${sample}_Aligned.sortedByCoord.out.bam"
    if [ ! -f "$bam" ]; then
        echo "[WARN] No se ha encontrado el BAM esperado: $bam. Se salta esta muestra."
        continue
    fi
    echo "    BAM: $bam"

    out_gtf="$SECOND_DIR/${sample}.gtf"

    stringtie "$bam" \
      -G "$MERGED_GTF" \
      -e -B -p "$THREADS" \
      -o "$out_gtf"

    # Añadimos al sample_list.csv para prepDE.py
    echo "${sample},${out_gtf}" >> sample_list.csv
done

##############################
# 4) MATRICES DE CONTEOS CON prepDE.py
##############################

echo "[INFO] Ejecutando prepDE.py..."
python "$PREPDE_SCRIPT" -i sample_list.csv

echo "[OK] Generados gene_count_matrix.csv y transcript_count_matrix.csv"

##############################
# 5) SOLAPAMIENTO TRANSCRITOS–TE USANDO SOLO GTF
##############################

# 5.1. Exones de merged.gtf (transcritos observados)
echo "[INFO] Extrayendo exones de merged.gtf..."
awk '$3=="exon"' "$MERGED_GTF" > merged_exons.gtf

# 5.2. Exones de repeatmasker.gtf (features de TE)
echo "[INFO] Extrayendo exones de repeatmasker.gtf..."
awk '$3=="exon"' "$RM_GTF" > repeat_exons.gtf

# 5.3. bedtools intersect directamente sobre GTF (-wa -wb para conservar atributos)
echo "[INFO] Intersectando exones de transcritos con exones de TE..."
bedtools intersect -wa -wb \
  -a merged_exons.gtf \
  -b repeat_exons.gtf \
  > TE_exon_overlap.gtf

echo "[OK] Generado TE_exon_overlap.gtf (GTF+GTF con info de solapamiento)"

cat << 'EOF'
[RESUMEN]
- merged.gtf                 -> catálogo global de transcritos
- gffcmp.annotated.gtf       -> class_code (novel vs genes)
- ${SECOND_DIR}/*.gtf        -> 2ª pasada de StringTie (un GTF por muestra)
- gene_count_matrix.csv
- transcript_count_matrix.csv
- merged_exons.gtf           -> exones de transcritos (StringTie)
- repeat_exons.gtf           -> exones de TE
- TE_exon_overlap.gtf        -> lineas exon_tx + exon_TE solapados
EOF