#!/bin/bash

sample_id="$1"
trimmed1="$2"
trimmed2="$3"
threads="$4"
MAP_DIR="$5"
GENOME_DIR="$6"
mismatch_align="$7"

# Setting up colors for messages

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/logging.sh"


# Script to process multiples SRA files creating results folders

 
    msg_info "==============================================="
    msg_info "Processing sample $SAMPLE_COUNT: $sample_id"
    msg_info "==============================================="
    
   
    # Verify fastq files do exist
    if [ -f "trim/${sample_id}_trimmed_1.fastq" ] && [ -f "trim/${sample_id}_trimmed_2.fastq" ]; then
        trimmed1="trim/${sample_id}_trimmed_1.fastq"
        trimmed2="trim/${sample_id}_trimmed_2.fastq"
    fi
    
    msg_info "[STAR] Using files: $trimmed1 and $trimmed2"
    
    # Verify genome index does exist
    if [ ! -d "../../$genome_index" ]; then
        msg_error "[STAR] The directory of genome index does not exist: genome_index"
        exit 1
    fi
    
    # STAR alignment
    msg_info "[STAR] Performing alignment for $sample_id..."
    STAR --runThreadN $threads \
      --genomeDir "$GENOME_DIR" \
      --readFilesIn "$trimmed1" "$trimmed2" \
      --outFileNamePrefix "$MAP_DIR/${sample_id}_" \
      --outSAMtype BAM SortedByCoordinate \
      --outFilterMultimapNmax 10 \
      --winAnchorMultimapNmax 50 \
      --outFilterMismatchNoverLmax $mismatch_align \
      --outSAMmultNmax 10 \
      --outMultimapperOrder Random \
      --runRNGseed 42 \
      --outSAMattributes NH HI AS nM NM MD

    
    # Verify BAM was created
    if [ ! -f "$MAP_DIR/${sample_id}_Aligned.sortedByCoord.out.bam" ]; then
        msg_error "[STAR] Alignment did not generate the expected BAM"
        cd "$CWD"
    fi

    samtools flagstat $MAP_DIR/${sample_id}_Aligned.sortedByCoord.out.bam
    msg_ok "[STAR] $SAMPLE_COUNT: $sample_id alignment completed"


