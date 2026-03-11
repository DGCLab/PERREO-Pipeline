#!/usr/bin/env bash
set -euo pipefail

CWD_OUTPUT="$(pwd)"
CWD="$(dirname "$CWD_OUTPUT")"

SAMPLESHEET="${CWD_OUTPUT}/samplesheet3.txt"
SAMPLES_DIR="${CWD_OUTPUT}/samples"

TX_FA="${CWD_OUTPUT}/nanopore_transcripts_unique.fa"   # <-- AJUSTA
SALMON_INDEX="${CWD_OUTPUT}/salmon_index"

THREADS=32

# 1) Índice
if [ -d "$SALMON_INDEX" ]; then
  echo "[SALMON] Index already exists in ${SALMON_INDEX}, creation omitted."
else
  echo "[SALMON] Building index from ${TX_FA}..."
  mkdir -p "$SALMON_INDEX"
  salmon index -t "$TX_FA" -i "$SALMON_INDEX" -k 31 --keepDuplicates

fi

infer_libtype () {
  local s="$1"
  s=$(echo "$s" | tr '[:upper:]' '[:lower:]' | xargs)
  case "$s" in
    forward)    echo "ISF" ;;
    reverse)    echo "ISR" ;;
    unstranded|none|na|"") echo "A" ;;
    *)          echo "A" ;;
  esac
}

echo "[SALMON] Starting quantification for all samples"

awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1, $2, $3}' "$SAMPLESHEET" \
  | while IFS=$'\t' read -r SAMPLE_ID STRANDEDNESS CONDITION; do
      [[ -z "$SAMPLE_ID" ]] && continue

      SAMPLE_DIR="${SAMPLES_DIR}/${SAMPLE_ID}"
      TRIM_DIR="${SAMPLE_DIR}/trim"

      R1=$(ls "${TRIM_DIR}/${SAMPLE_ID}_trimmed_1.fastq"* 2>/dev/null || true)
      R2=$(ls "${TRIM_DIR}/${SAMPLE_ID}_trimmed_2.fastq"* 2>/dev/null || true)

      if [[ -z "$R1" || -z "$R2" ]]; then
        echo "[SALMON] WARNING: trimmed FASTQ not found for ${SAMPLE_ID} in ${TRIM_DIR}, skipping." >&2
        continue
      fi

      SAMPLE_OUT="${SAMPLE_DIR}/salmon_results"
      if [[ -d "$SAMPLE_OUT" && -f "${SAMPLE_OUT}/quant.sf" ]]; then
        echo "[SALMON] Skipping ${SAMPLE_ID}, quantification already exists."
        continue
      fi

      LIBTYPE=$(infer_libtype "$STRANDEDNESS")
      echo "[SALMON] Quantifying ${SAMPLE_ID} (strandedness=${STRANDEDNESS} -> libType=${LIBTYPE})"

      mkdir -p "$SAMPLE_OUT"

      /usr/bin/time -v -o "${SAMPLE_OUT}/${SAMPLE_ID}_salmon_time.log" \
      salmon quant \
        -i "$SALMON_INDEX" \
        -l "$LIBTYPE" \
        -1 "$R1" \
        -2 "$R2" \
        -p "$THREADS" \
        --validateMappings \
        --numBootstraps 30 \
        --seqBias \
        --gcBias \
        --posBias \
        -o "$SAMPLE_OUT"
      
    done

echo "[SALMON] Quantification for all samples completed."
