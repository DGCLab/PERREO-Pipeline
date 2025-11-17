#!/usr/bin/env bash
set -euo pipefail

#Defining current path
CWD="$(pwd)"
GLOBAL_DIR="$CWD"

# Scripts paths that already exist (ajusta nombres/paths)
TRIM_SCRIPT="$CWD/scripts_SE/script_trimming_SE.sh"   
MAP_SCRIPT="$CWD/scripts_SE/script_T2T_star_map_se.sh"
MARKDUP_SCRIPT="$CWD/scripts_SE/script_markduplicates_se.sh"           
QUANT_SCRIPT="$CWD/scripts_SE/quant.R"
MERGE_QUANT_SCRIPT="$CWD/scripts_SE/merge_quant.R"
DEA_SCRIPT="$CWD/scripts_SE/dea.R"      
PRED_MODEL="$CWD/scripts_SE/prediction_model.R"     
STAR_PATH="$CWD/star-2.7.10b-h9ee0642_0/bin/STAR"
PICARD_PATH="$CWD/picard.jar"

sample_list=""
reference_genome=""
genome_gtf=""
repeat_gtf=""
threads=8
adapter=""
mismatch_align=0.05
project_name=""
remove_duplicates=""
batch_effect=""
method=""
k_num=2


  # Parsing arguments
while [[ $# -gt 0 ]]; do
  case $1 in
      -sample_list) sample_list="$2"; shift 2 ;;
      -reference_genome) reference_genome="$2"; shift 2 ;;
      -genome_gtf) genome_gtf="$2"; shift 2 ;;
      -repeat_gtf) repeat_gtf="$2"; shift 2 ;;
      -threads) threads="$2"; shift 2 ;;
      -adapter) adapter="$2"; shift 2 ;;
      -mismatch_align) mismatch_align="$2"; shift 2 ;;
      -project_name) project_name="$2"; shift 2 ;;
      -remove_duplicates) remove_duplicates="$2"; shift 2 ;;
      -batch_effect) batch_effect="$2"; shift 2 ;;
      -method) method="$2"; shift 2 ;;
      *) echo "Opción desconocida: $1"; shift ;;
  esac
done


# Function for running the analysis with all the samples

run_pipeline_sample() {
  echo "samplesheet: $sample_list"
  echo "Reference: $reference_genome"
  echo "GTF: $genome_gtf"
  echo "Repeats: $repeat_gtf"
  echo "Threads: $threads"
  echo "adaptor: $adapter"
  echo "mismatch_align: $mismatch_align"
  echo "Proyect: $project_name"
  echo "Remove duplicates: $remove_duplicates"
  echo "Batch effect: $batch_effect"
  echo "Method: $method"

if [ -d "genome_index" ]; then
  echo "✅ Genome index already exists in genome_index, creation omitted."
else
  echo "Creating genome index..."
  "$STAR_PATH" --runMode genomeGenerate --runThreadN $threads --genomeDir genome_index --genomeFastaFiles "$reference_genome" --sjdbGTFfile "$genome_gtf" --sjdbOverhang 100
fi

GENOME_DIR=$(realpath genome_index)

cd SAMPLES
SAMPLES_DIR=$CWD/SAMPLES


# Requiere: CWD y sample_list definidos, y TRIM_* scripts disponibles
set -euo pipefail

awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1, $2, $3}' "../$sample_list" \
| while IFS=$'\t' read -r sample_id STRAND CONDITION; do
    [[ -z "$sample_id" ]] && continue

    cutadapt_threads=$threads

    STRAND=$(echo "$STRAND" | tr '[:upper:]' '[:lower:]' | xargs)
    CONDITION=$(echo "${CONDITION:-}" | xargs)

    echo "Processing $sample_id with strandedness=$STRAND and condition=$CONDITION"

    # Rutas por muestra
    SAMPLE_DIR="$CWD/SAMPLES/${sample_id}"
    TRIM_DIR="${SAMPLE_DIR}/trim"
    MAP_DIR="${SAMPLE_DIR}/alignment"
    IN="${sample_id}.fastq"

    # Asegura estructura básica
    mkdir -p "$SAMPLE_DIR"

    # Entra en la carpeta de la muestra
    cd "$SAMPLE_DIR" || { echo "Error: Unable to enter in $SAMPLE_DIR"; exit 1; }

    # Comprueba inputs crudos
    if [ ! -f "$IN" ]; then
        echo "⚠️  Raw files not found for $sample_id (IN=$IN)"
        continue
    fi
    echo "IN: $IN"

    echo "STRAND: [$STRAND]"

    # Sólo hace trimming si faltan los FASTQ trimmados
    if [ ! -f "$TRIM_DIR/${sample_id}_trimmed.fastq" ]; then
        mkdir -p "$TRIM_DIR" "$MAP_DIR"

        # ------------ 1) TRIMMING READS ---------------------------
        bash "$TRIM_SCRIPT" "$sample_id" "$IN" "$TRIM_DIR" "$adapter" "$threads"
    else
        echo "✓ Trimmed FASTQs already exist for $sample_id — skipping trimming."
    fi

done


cd $SAMPLES_DIR


# ---------- 2) MAPPING AGAINST REFERENCE GENOME ----------

  awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1, $2, $3}' "../$sample_list" \
  | while IFS=$'\t' read -r sample_id STRAND CONDITION; do
     [[ -z "$sample_id" ]] && continue

     STRAND=$(echo "$STRAND" | tr '[:upper:]' '[:lower:]' | xargs)
     SAMPLE_DIR="$CWD/SAMPLES/${sample_id}"
     TRIM_DIR="${SAMPLE_DIR}/trim"
     MAP_DIR="${SAMPLE_DIR}/alignment"

  if [[ -f "$MAP_DIR/${sample_id}_Aligned.sortedByCoord.out.bam" ]];then
        echo 'Skipping mapping...'
  else
        bash "$MAP_SCRIPT" "$sample_id" "${TRIM_DIR}/${sample_id}_trimmed.fastq" "$threads" "$MAP_DIR" "$STAR_PATH" "$GENOME_DIR" "$mismatch_align"
  fi
  done


# -----------3) MARKDUPLICATES ----------------------------

  
 awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1, $2, $3}' "../$sample_list" \
  | while IFS=$'\t' read -r sample_id STRAND CONDITION; do
     [[ -z "$sample_id" ]] && continue
  
  SAMPLE_DIR="$CWD/SAMPLES/${sample_id}"
  MAP_DIR="${SAMPLE_DIR}/alignment"

  if [[ -f "$MAP_DIR/${sample_id}_marked_duplicates_STAR.bam" ]];then
        echo 'Skipping Markduplicates...'
  else
        REMOVE_DUPLICATES=false
        bash "$MARKDUP_SCRIPT" "$sample_id" "$threads" "$MAP_DIR" "$PICARD_PATH" "$remove_duplicates"
  fi
  done

#multiqc .


# ---------- 4) QUANTIFICATION ----------------------------

  awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1, $2, $3}' "../$sample_list" \
  | while IFS=$'\t' read -r sample_id STRAND CONDITION; do
      [[ -z "$sample_id" ]] && continue

      SAMPLE_DIR="$CWD/SAMPLES/${sample_id}"
      MAP_DIR="${SAMPLE_DIR}/alignment"

   if [ ! -d "${SAMPLE_DIR}/Quantification" ]; then 
      mkdir ${SAMPLE_DIR}/Quantification
   fi

      QUANT_DIR="${SAMPLE_DIR}/Quantification"

      STRAND=$(echo "$STRAND" | tr '[:upper:]' '[:lower:]' | xargs)

      REP_GTF_PATH=$CWD/$repeat_gtf 

   if [ ! -f "$QUANT_DIR/${sample_id}_quant.txt" ]; then

   #Running Rsubread for quantification
      Rscript "$QUANT_SCRIPT" "$MAP_DIR" "$sample_id" "$REP_GTF_PATH" "$threads" "$STRAND" "$SAMPLE_DIR" "$QUANT_DIR"
   fi

  done



# ---------- 5) COUNT MATRIXES MERGE ----------------------


#Finally, we call the last script to merge all the count matrixes

   if [ ! -d "$CWD/SAMPLES/DEA_results" ]; then

      mkdir $CWD/SAMPLES/DEA_results
   fi

   DEA_results="$CWD/SAMPLES/DEA_results"
      
   REP_GTF_PATH=$CWD/$repeat_gtf 

   if [ ! -f "$CWD/SAMPLES/count_data.txt" ]; then
     echo 'count_data.txt must be generated'

     Rscript "$MERGE_QUANT_SCRIPT" "$GLOBAL_DIR" "$SAMPLES_DIR" "$REP_GTF_PATH" "$threads" "$DEA_results"
     cd ..
   fi



# ---------- 6) DIFFERENTIAL EXPRESSION ANALYSIS ----------

Rscript "$DEA_SCRIPT" "$batch_effect" "$sample_list" "$method" "$CWD" "$REP_GTF_PATH" "$k_num"


# ---------- 7) PREDICTION MODEL ANALYSIS -----------------

rows=$(( $(wc -l < "$CWD/$sample_list") - 1 ))
if [ "$rows" -gt 5 ]; then
Rscript "$PRED_MODEL" "$CWD" "$sample_list"
fi

}


#Running the function

run_pipeline_sample 
