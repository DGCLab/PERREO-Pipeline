#!/bin/bash

sample_id="$1"
threads="$2"
MAP_DIR="$3"
remove_duplicates="$4"


       # Add read group with AddOrReplaceReadGroups
       echo "Adding read group for $sample_id..."
       picard AddOrReplaceReadGroups \
        I="$MAP_DIR/${sample_id}_Aligned.sortedByCoord.out.bam" \
        O="$MAP_DIR/${sample_id}_with_readgroup_STAR.bam" \
        RGID="group${SAMPLE_COUNT}" \
        RGLB=lib1 \
        RGPL=illumina \
        RGPU=unit1 \
        RGSM="sample${SAMPLE_COUNT}"
    
    
      # Marking duplicates 
      echo "Marcando duplicados para $sample_id..."
       picard MarkDuplicates \
         I="$MAP_DIR/${sample_id}_with_readgroup_STAR.bam" \
         O="$MAP_DIR/${sample_id}_marked_duplicates_STAR.bam" \
         M="$MAP_DIR/${sample_id}_marked_dup_metrics_STAR.txt" \
         REMOVE_DUPLICATES=$remove_duplicates

       echo "Procesamiento completado para $sample_id"
    
       # Returning to base directory for the next iteration
       cd "$CWD"

    
