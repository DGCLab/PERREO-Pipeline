#!/bin/bash

sample_id="$1"
threads="$2"
MAP_DIR="$3"
PICARD_PATH="$4"
REMOVE_DUPLICATES="$5"


       # Añadir grupo de lectura con Picard
       echo "Añadiendo grupo de lectura con Picard para $sample_id..."
       java -jar "$PICARD_PATH" AddOrReplaceReadGroups \
        I="$MAP_DIR/${sample_id}_Aligned.sortedByCoord.out.bam" \
        O="$MAP_DIR/${sample_id}_with_readgroup_STAR.bam" \
        RGID="group${SAMPLE_COUNT}" \
        RGLB=lib1 \
        RGPL=illumina \
        RGPU=unit1 \
        RGSM="sample${SAMPLE_COUNT}"
    
    
      # Eliminación de duplicados
      echo "Eliminando duplicados para $sample_id..."
       java -jar "$PICARD_PATH" MarkDuplicates \
         I="$MAP_DIR/${sample_id}_with_readgroup_STAR.bam" \
         O="$MAP_DIR/${sample_id}_marked_duplicates_STAR.bam" \
         M="$MAP_DIR/${sample_id}_marked_dup_metrics_STAR.txt" \
         REMOVE_DUPLICATES=$REMOVE_DUPLICATES
    
       echo "Procesamiento completado para $sample_id"
    
       # Volver al directorio base para la siguiente iteración
       cd "$CWD"
    
       # Incrementar el contador de muestras
       #SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
    
