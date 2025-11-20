#!/bin/bash

sample_id="$1"
trimmed1="$2"
trimmed2="$3"
threads="$4"
MAP_DIR="$5"
GENOME_DIR="$6"
mismatch_align="$7"

# Script to process multiples SRA files creating results folders

 
    echo "==============================================="
    echo "Procesando muestra $SAMPLE_COUNT: $sample_id"
    echo "==============================================="
    
   
    # Verify fastq files do exist
    if [ -f "trim/${sample_id}_trimmed_1.fastq" ] && [ -f "trim/${sample_id}_trimmed_2.fastq" ]; then
        trimmed1="trim/${sample_id}_trimmed_1.fastq"
        trimmed2="trim/${sample_id}_trimmed_2.fastq"
    fi
    
    echo "Usando archivos: $trimmed1 y $trimmed2"
    
    # Verify genome index does exist
    if [ ! -d "../../$genome_index" ]; then
        echo "ERROR: The directory of genome index does not exist: genome_index"
        exit 1
    fi
    
    # STAR alignment
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

    
    # Verify BAM was created
    if [ ! -f "$MAP_DIR/${sample_id}_Aligned.sortedByCoord.out.bam" ]; then
        echo "ERROR: STAR alignment did not generate the expected BAM"
        cd "$CWD"
    fi

    samtools flagstat $MAP_DIR/${sample_id}_Aligned.sortedByCoord.out.bam


echo "Alignment for all samples completed"
