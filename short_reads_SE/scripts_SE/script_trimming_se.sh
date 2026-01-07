#!/bin/bash

sample_id="$1"
IN="$2"
TRIM_DIR="$3"
adapter="$4"
trimming_type="$5"
threads="$6"
trimming_quality="$7"
min_length="$8"
initial_trim_read="$9"
polya="${10}"

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
if [ -f "trim/${sample_id}_trimmed.fastq" ]; then
    msg_ok "[CUTADAPT] Files already exist in trim folder, creation omitted."
else
  if [ -z "$adapter" ]; then
    msg_info '[CUTADAPT] Performing trimming without removing adapters (already removed)'
  cutadapt -j "$threads" -q "$trimming_quality","$trimming_quality" \
       --trim-n -m "$min_length" -u "$initial_trim_read" ${polya:+$polya} \
        -o "${TRIM_DIR}/${sample_id}_trimmed.fastq" \
        "$IN"  > cutadapt.log 2>&1
  else
     msg_info '[CUTADAPT] Performing trimming removing adapters'
    cutadapt -j "$threads" -q "$trimming_quality","$trimming_quality" \
        -a "$adapter" --trim-n -m "$min_length" -u "$initial_trim_read" ${polya:+$polya} \
        -o "${TRIM_DIR}/${sample_id}_trimmed.fastq" \
        "$IN"  > cutadapt.log 2>&1
fi

fastqc ${TRIM_DIR}/${sample_id}_trimmed.fastq

