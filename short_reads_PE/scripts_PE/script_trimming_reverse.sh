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

echo $sample_id
echo $IN1
echo $IN2
echo $TRIM_DIR
echo $adapt_r1
echo $adapt_r2
echo $trimming_type
echo $threads
echo $trimming_quality
echo $min_length

#Run cutadapt
if [[ -f "trim/${sample_id}_trimmed_1.fastq.gz" && -f "trim/${sample_id}_trimmed_2.fastq.gz" ]]; then
    echo "✅ Files already exist in trim folder, creation omitted."
else
  if [[ -z "$adapt_r1" && -z "$adapt_r2" ]]; then
    echo 'Performing trimming without removing adapters(already removed)'
    cutadapt -j "$threads" --pair-filter any -q "$trimming_quality","$trimming_quality" \
        --trim-n -m "$min_length" \
        -o "${TRIM_DIR}/${sample_id}_trimmed_1.fastq" \
        -p "${TRIM_DIR}/${sample_id}_trimmed_2.fastq" \
        "$IN1" "$IN2" > cutadapt.log 2>&1
  else
    echo 'Performing trimming removing adapters'
   cutadapt -j "$threads" --pair-filter any -q "$trimming_quality","$trimming_quality" \
        -a "$adapt_r1" -A "$adapt_r2" --trim-n -m "$min_length" \
        -o "${TRIM_DIR}/${sample_id}_trimmed_1.fastq" \
        -p "${TRIM_DIR}/${sample_id}_trimmed_2.fastq" \
        "$IN1" "$IN2" > cutadapt.log 2>&1

  fi



  if [[ "$trimming_type" == "trimming_extra" ]]; then
       echo "→ Performing extra trimming"
       gzip "${TRIM_DIR}/${sample_id}_trimmed_1.fastq"
       gzip "${TRIM_DIR}/${sample_id}_trimmed_2.fastq"
           
       # Verifying that the fastq files already exists
       if [ -f "trim/${sample_id}_trimmed_1.fastq.gz" ] && [ -f "trim/${sample_id}_trimmed_2.fastq.gz" ]; then
          python $TRIM_EXTRA -s reverse -o "trim/${sample_id}_trimmed" -i "trim/${sample_id}_trimmed" > trimGC.log 2>&1
          gunzip ${TRIM_DIR}/${sample_id}_trimmed_1.fastq.gz
          gunzip ${TRIM_DIR}/${sample_id}_trimmed_2.fastq.gz


       else
          echo "ERROR: Trimmed files not found"
       fi 
   
  else
       echo "→ Skipping extra trimming"
       
  fi
fi

fastqc ${TRIM_DIR}/${sample_id}_trimmed_1.fastq
fastqc ${TRIM_DIR}/${sample_id}_trimmed_2.fastq
