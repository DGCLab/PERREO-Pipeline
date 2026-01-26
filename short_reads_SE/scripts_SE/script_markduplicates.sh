#!/bin/bash

sample_id="$1"
threads="$2"
MAP_DIR="$3"
remove_duplicates="$4"

# Setting up colors for messages

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/logging.sh"

###

       # Add read group with AddOrReplaceReadGroups
       msg_info "[MARKDUPLICATES] Adding read group for $sample_id..."
       picard AddOrReplaceReadGroups \
        I="$MAP_DIR/${sample_id}_Aligned.sortedByCoord.out.bam" \
        O="$MAP_DIR/${sample_id}_with_readgroup_STAR.bam" \
        RGID="group${SAMPLE_COUNT}" \
        RGLB=lib1 \
        RGPL=illumina \
        RGPU=unit1 \
        RGSM="sample${SAMPLE_COUNT}"
    
    
      # Marking duplicates 
      msg_info "[MARKDUPLICATES] Marking duplicates for $sample_id..."
       picard -Xmx"${ram}" MarkDuplicates \
         I="$MAP_DIR/${sample_id}_with_readgroup_STAR.bam" \
         O="$MAP_DIR/${sample_id}_marked_duplicates_STAR.bam" \
         M="$MAP_DIR/${sample_id}_marked_dup_metrics_STAR.txt" \
         REMOVE_DUPLICATES=$remove_duplicates

       samtools sort -@ 8 $MAP_DIR/${sample_id}_marked_duplicates_STAR.bam -o $MAP_DIR/${sample_id}_marked_duplicates_STAR_sorted.bam
       samtools index $MAP_DIR/${sample_id}_marked_duplicates_STAR_sorted.bam

       msg_ok "[MARKDUPLICATES] Processing completed for $sample_id"
    
       # Returning to base directory for the next iteration
       cd "$CWD"

    
