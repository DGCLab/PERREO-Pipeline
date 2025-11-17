#!/bin/bash

sample_id="$1"
trimmed1="$2"
trimmed2="$3"
threads="$4"
MAP_DIR="$5"
#STAR_PATH="$6"
GENOME_DIR="$6"
mismatch_align="$7"

# Script para procesamiento de multiples archivos SRA con creación de carpetas de resultados
 
    echo "==============================================="
    echo "Procesando muestra $SAMPLE_COUNT: $sample_id"
    echo "==============================================="
    
   
    # Verificar que los archivos FASTQ existen
    if [ -f "trim/${sample_id}_trimmed_1.fastq" ] && [ -f "trim/${sample_id}_trimmed_2.fastq" ]; then
        trimmed1="trim/${sample_id}_trimmed_1.fastq"
        trimmed2="trim/${sample_id}_trimmed_2.fastq"
    #else
      # gunzip trim/${sample_id}_trimmed_1.fastq.gz
       # gunzip trim/${sample_id}_trimmed_2.fastq.gz
        #trimmed1="trim/${sample_id}_trimmed_1.fastq"
        #trimmed2="trim/${sample_id}_trimmed_2.fastq"
    fi
    
    echo "Usando archivos: $trimmed1 y $trimmed2"
    
    # Verificar que el índice del genoma existe
    if [ ! -d "../../$genome_index" ]; then
        echo "ERROR: El directorio del índice del genoma no existe: genome_index"
        exit 1
    fi
    
    # Alineamiento con STAR - los resultados se guardan en la carpeta específica de sample_id dentro de Nanopore
    echo "Realizando alineamiento con STAR para $sample_id..."
    STAR --runThreadN $threads \
      --genomeDir "$GENOME_DIR" \
      --readFilesIn "$trimmed1" "$trimmed2" \
      --outFileNamePrefix "$MAP_DIR/${sample_id}_" \
      --outSAMtype BAM SortedByCoordinate \
      --outFilterMultimapNmax 500 \
      --winAnchorMultimapNmax 500 \
      --outFilterMismatchNoverLmax $mismatch_align \
      --outSAMmultNmax 500 \
      --outMultimapperOrder Random \
      --runRNGseed 42 \
      --outSAMattributes NH HI AS nM NM MD

    
    # Verificar si se creó el BAM
    if [ ! -f "$MAP_DIR/${sample_id}_Aligned.sortedByCoord.out.bam" ]; then
        echo "ERROR: El alineamiento con STAR no generó el archivo BAM esperado"
        cd "$CWD"
        #continue
    fi

    samtools flagstat $MAP_DIR/${sample_id}_Aligned.sortedByCoord.out.bam


echo "Alineamiento de todas las muestras completado"