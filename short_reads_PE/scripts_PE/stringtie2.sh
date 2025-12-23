#!/bin/bash

# Setting up colors for messages

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/logging.sh"

combined_annotation="$1"
sample="$2"
threads="$3"
strandedness="$4"
CWD="$5"

# Checking samples folders
if [ $# -eq 0 ]; then
    msg_error "[STRINGTIE2] You didn’t specify any directories."
    exit 1
fi

# Running above the samples folders selected
    # Checking the folder exists

    # Buscar el BAM dentro de la carpeta
    BAM=$(find "$CWD/samples/$sample/alignment" -maxdepth 1 -name "*marked_duplicates_STAR.bam" | head -n 1)

    if [ -z "$BAM" ]; then
        msg_error "[STRINGTIE2] No .bam found inside $sample"
    fi

    msg_ok "[STRINGTIE2] BAM found: $BAM"

    # Basename of the folder
    BASENAME=$(basename "$sample")

    # Output name
    OUTPUT_GTF="$CWD/samples/$sample/${BASENAME}_transcriptome.gtf"

    msg_info "[STRINGTIE2] Running StringTie..."

    if [ "$strandedness" = "forward" ]; then
    msg_info "[STRINGTIE2] Using forward parameter"
    stringtie "$BAM" \
        -L \
        -p "$threads" \
        -G "$combined_annotation" \
        -o "$OUTPUT_GTF" \
        --fr \
	-l perreo
    fi

    if [ "$strandedness" = "reverse" ]; then
    msg_info "[STRINGTIE2] Using reverse parameter"
    stringtie "$BAM" \
        -L \
        -p "$threads" \
        -G "$combined_annotation" \
        -o "$OUTPUT_GTF" \
        --rf \
	-l perreo
    fi

    msg_ok "[STRINGTIE2] Transcriptome generated: $OUTPUT_GTF"

msg_info "---------------------------------------------"
msg_info "   FINISHED"
msg_info "---------------------------------------------"
