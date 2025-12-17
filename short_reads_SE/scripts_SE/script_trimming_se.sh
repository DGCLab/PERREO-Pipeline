#!/bin/bash

sample_id="$1"
IN="$2"
TRIM_DIR="$3"
adapter="$4"
trimming_type="$5"
threads="$6"
trimming_quality="$7"
min_length="$8"

# Setting up colors for messages

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/logging.sh"

####

msg_info $sample_id
msg_info $IN
msg_info $TRIM_DIR
msg_info $adapter
msg_info $trimming_type
msg_info $threads
msg_info $trimming_quality
msg_info $min_length

#Corro cutadapt
    # It is necessary to activate the conda environment where cutadapt software is installed
if [ -f "trim/${sample_id}_trimmed.fastq ]; then
    echo "✅ Archivos ya existen en trim/, se omite la creación."
else
    cutadapt -j "$threads" --pair-filter any -q "$trimming_quality","$trimming_quality" \
        -a "$adapter" --trim-n -m "$min_length" \
        -o "${TRIM_DIR}/${sample_id}_trimmed.fastq" \
        "$IN"  > cutadapt.log 2>&1
       
  fi
fi

fastqc ${TRIM_DIR}/${sample_id}_trimmed.fastq

