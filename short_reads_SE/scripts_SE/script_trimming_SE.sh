#!/bin/bash

sample_id="$1"
IN="$2"
TRIM_DIR="$3"
adapter="$4"
threads="$5"

echo $sample_id
echo $IN
echo $TRIM_DIR
echo $adapter
echo $threads

#Corro cutadapt
    # It is necessary to activate the conda environment where cutadapt software is installed
if [ -f "trim/${sample_id}_trimmed.fastq ]; then
    echo "✅ Archivos ya existen en trim/, se omite la creación."
else
    cutadapt -j "$threads" --pair-filter any -q 30,30 \
        -a "$adapter" --trim-n -m 16 \
        -o "${TRIM_DIR}/${sample_id}_trimmed.fastq" \
        "$IN"  > cutadapt.log 2>&1
       
  fi
fi

fastqc ${TRIM_DIR}/${sample_id}_trimmed.fastq

