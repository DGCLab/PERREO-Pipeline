
#!/bin/bash

sample_id="$1"
IN1="$2"
IN2="$3"
TRIM_DIR="$4"
adapt_r1="$5"
adapt_r2="$6"
trimming_type="$7"
threads="$8"
trimming_quality="$9"
min_length="${10}"
initial_trim_read1="${11}"
initial_trim_read2="${12}"
polya="${13}"
TRIM_EXTRA="${14}"


# Setting up colors for messages

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/logging.sh"

####

msg_info $sample_id
msg_info $IN1
msg_info $IN2
msg_info $TRIM_DIR
msg_info $adapt_r1
msg_info $adapt_r2
msg_info $trimming_type
msg_info $threads
msg_info $trimming_quality
msg_info $min_length
msg_info $initial_trim_read1
msg_info $initial_trim_read2

#Run cutadapt
if [[ -f "trim/${sample_id}_trimmed_1.fastq.gz" && -f "trim/${sample_id}_trimmed_2.fastq.gz" ]]; then
    msg_ok "[CUTADAPT] Files already exist in trim folder, creation omitted."
else
  if [[ -z "$adapt_r1" && -z "$adapt_r2" ]]; then
    msg_info '[CUTADAPT] Performing trimming without removing adapters (already removed)'
    cutadapt -j "$threads" --pair-filter any -q "$trimming_quality","$trimming_quality" \
        --trim-n -m "$min_length" -u "$initial_trim_read1"  -U "$initial_trim_read2" ${polya:+$polya}  \
        -o "${TRIM_DIR}/${sample_id}_trimmed_1.fastq" \
        -p "${TRIM_DIR}/${sample_id}_trimmed_2.fastq" \
        "$IN1" "$IN2" > cutadapt.log 2>&1
  else
   msg_info '[CUTADAPT] Performing trimming removing adapters'
   cutadapt -j "$threads" --pair-filter any -q "$trimming_quality","$trimming_quality" \
        -a "$adapt_r1" -A "$adapt_r2" --trim-n -m "$min_length" -u "$initial_trim_read1"  -U "$initial_trim_read2" ${polya:+$polya}  \
        -o "${TRIM_DIR}/${sample_id}_trimmed_1.fastq" \
        -p "${TRIM_DIR}/${sample_id}_trimmed_2.fastq" \
        "$IN1" "$IN2" > cutadapt.log 2>&1

  fi



  if [[ "$trimming_type" == "extra" ]]; then
       msg_info "[CUTADAPT] → Performing extra trimming"
       gzip "${TRIM_DIR}/${sample_id}_trimmed_1.fastq"
       gzip "${TRIM_DIR}/${sample_id}_trimmed_2.fastq"
           
       # Verifying that the fastq files already exists
       if [ -f "trim/${sample_id}_trimmed_1.fastq.gz" ] && [ -f "trim/${sample_id}_trimmed_2.fastq.gz" ]; then
          python "$TRIM_EXTRA" -s reverse -o "trim/${sample_id}_trimmed_gc" -i "trim/${sample_id}_trimmed" > trim_extra.log 2>&1
          gunzip ${TRIM_DIR}/${sample_id}_trimmed_gc_1.fastq.gz
          gunzip ${TRIM_DIR}/${sample_id}_trimmed_gc_2.fastq.gz
          
          mv ${TRIM_DIR}/${sample_id}_trimmed_gc_1.fastq ${TRIM_DIR}/${sample_id}_trimmed_1.fastq
          mv ${TRIM_DIR}/${sample_id}_trimmed_gc_2.fastq ${TRIM_DIR}/${sample_id}_trimmed_2.fastq


       else
          msg_error "[CUTADAPT] Trimmed files not found"
       fi 
   
  else
       msg_info "[CUTADAPT] → Skipping extra trimming"
       
  fi
fi

fastqc ${TRIM_DIR}/${sample_id}_trimmed_1.fastq
fastqc ${TRIM_DIR}/${sample_id}_trimmed_2.fastq
