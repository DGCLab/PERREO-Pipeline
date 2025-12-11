#!/bin/bash


# ANOTACIÓN (modifica si quieres otra)
CWD="$1"
combined_annotation="$2"
sample="$3"
threads="$4"
strandedness="$5"

echo "---------------------------------------------"
echo "   Running StringTie2 Long-Reads Assembly"
echo "---------------------------------------------"
echo "Annotation: $combined_annotation"
echo "Threads: $threads"
echo ""

# Checking samples folders
if [ $# -eq 0 ]; then
    echo "ERROR: No has especificado carpetas."
    echo "Uso: bash run_stringtie.sh carpeta1 carpeta2 ..."
    exit 1
fi

# Running above the samples folders selected
    # Checking the folder exists

    # Buscar el BAM dentro de la carpeta
    BAM=$(find "$sample" -maxdepth 1 -name "*.bam" | head -n 1)

    if [ -z "$BAM" ]; then
        echo "  ERROR: NO .bam found inside $sample"
    fi

    echo "  BAM found: $BAM"

    # Basename of the folder
    BASENAME=$(basename "$sample")

    # Output name
    OUTPUT_GTF="$sample/${BASENAME}_transcriptome.gtf"

    echo "  Running StringTie..."
     if [ "$strandedness" = "forward" ]; then
    echo "Using forward parameter"
    stringtie "$BAM" \
        -L \
        -p "$threads" \
        -G "$combined_annotation" \
        -o "$OUTPUT_GTF" \
        --fr \
	-l perreo
    fi

    if [ "$strandedness" = "reverse" ]; then
    echo "Using reverse parameter"
    stringtie "$BAM" \
        -L \
        -p "$threads" \
        -G "$combined_annotation" \
        -o "$OUTPUT_GTF" \
        --rf \
	-l perreo
    fi

    echo "  ✔ Transcriptome generated: $OUTPUT_GTF"
    echo ""

echo "---------------------------------------------"
echo "   FINISHED"
echo "---------------------------------------------"
