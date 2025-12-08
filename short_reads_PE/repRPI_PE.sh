#!/usr/bin/env bash
set -euo pipefail

#Defining current path
CWD="$(pwd)"

# Scripts paths that already exist (ajusta nombres/paths)
TRIM_FW_SCRIPT="$CWD/scripts_PE/script_trimming_forward.sh"   
TRIM_RV_SCRIPT="$CWD/scripts_PE/script_trimming_reverse.sh"
TRIM_EXTRA="$CWD/scripts_PE/trimGC.py"
MAP_SCRIPT="$CWD/scripts_PE/script_T2T_star_map.sh"
MARKDUP_SCRIPT="$CWD/scripts_PE/script_markduplicates.sh"           
QUANT_SCRIPT="$CWD/scripts_PE/quant.R"
MERGE_QUANT_SCRIPT="$CWD/scripts_PE/merge_quant.R"
DEA_SCRIPT="$CWD/scripts_PE/dea.R"
DEA_SCRIPT_multicond="$CWD/scripts_PE/dea_multicond.R"
ASSEMBLY_SCRIPT="$CWD/scripts_PE/stringtie2.sh"   
ASSEMBLY_SCRIPT_2="$CWD/scripts_PE/stringtie2_2.0.sh" 
PREPDE_SCRIPT="$CWD/scripts_PE/prepDE.py3" 
WGCNA_SCRIPT="$CWD/scripts_PE/WGCNA.R"         
PRED_MODEL="$CWD/scripts_PE/prediction_model.R"     

  # Parsing arguments
while [[ $# -gt 0 ]]; do
  case $1 in
      -sample_list) sample_list="$2"; shift 2 ;;
      -reference_genome) reference_genome="$2"; shift 2 ;;
      -genome_gtf) genome_gtf="$2"; shift 2 ;;
      -repeat_gtf) repeat_gtf="$2"; shift 2 ;;
      -threads) threads="$2"; shift 2 ;;
      -adapt_r1) adapt_r1="$2"; shift 2 ;;
      -adapt_r2) adapt_r2="$2"; shift 2 ;;
      -trimming_quality_threshold) trimming_quality_threshold="$2"; shift 2 ;;
      -min_length_trim) min_length_trim="$2"; shift 2 ;;
      -max_length_trim) min_length_trim="$2"; shift 2 ;;
      -trimming_type) trimming_type="$2"; shift 2 ;;
      -mismatch_align) mismatch_align="$2"; shift 2 ;;
      -project_name) project_name="$2"; shift 2 ;;
      -remove_duplicates) remove_duplicates="$2"; shift 2 ;;
      -log2FC) log2FC="$2"; shift 2 ;;
      -FDR) FDR="$2"; shift 2 ;;
      -batch_effect) batch_effect="$2"; shift 2 ;;
      -method) method="$2"; shift 2 ;;
      -prediction_model) method="$2"; shift 2 ;;
      *) echo "Opción desconocida: $1"; shift ;;
  esac
done

# --- Assigning values by default ---
threads=${threads:-8}
k_num=${k_num:-2}
FDR=${FDR:-0.05}
log2FC=${log2FC:-1.0}
mismatch_align=${mismatch_align:-0.05}
trimming_quality_threshold=${trimming_quality_threshold:-30}
min_length_trim=${min_length_trim:-16}
prediction_model=${prediction_model:-yes}


# Function for running the analysis with all the samples

run_pipeline_sample() {
  echo " ----- Pipeline paramters summary -----"
  echo "samplesheet: $sample_list"
  echo "Reference: $reference_genome"
  echo "GTF: $genome_gtf"
  echo "Repeats: $repeat_gtf"
  echo "Threads: $threads"
  echo "R1 adaptor: $adapt_r1"
  echo "R2 adaptor: $adapt_r2"
  echo "trimming type: $trimming_type"
  echo "trimming_quality_threshold: $trimming_quality_threshold" 
  echo "minimum_length_trim: $min_length_trim"
  echo "mismatch_align type: $mismatch_align"
  echo "Proyect: $project_name"
  echo "Remove duplicates: $remove_duplicates"
  echo "Batch effect: $batch_effect"
  echo "Method: $method"
  echo "log2FC: $log2FC"
  echo "FDR: $FDR"

if [ -d "genome_index" ]; then
  echo "✅ Genome index already exists in genome_index, creation omitted."
else
  echo "Creating genome index..."
  STAR --runMode genomeGenerate --runThreadN $threads --genomeDir genome_index --genomeFastaFiles "$reference_genome" --sjdbGTFfile "$genome_gtf" --sjdbOverhang 100
fi

GENOME_DIR=$(realpath genome_index)

cd SAMPLES
SAMPLES_DIR=$CWD/SAMPLES

set -euo pipefail

awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1, $2, $3}' "../$sample_list" \
| while IFS=$'\t' read -r sample_id STRAND CONDITION; do
    [[ -z "$sample_id" ]] && continue

    cutadapt_threads=$threads

    STRAND=$(echo "$STRAND" | tr '[:upper:]' '[:lower:]' | xargs)
    CONDITION=$(echo "${CONDITION:-}" | xargs)

    echo "Processing $sample_id with strandedness=$STRAND and condition=$CONDITION"

    # Path per sample
    SAMPLE_DIR="$CWD/SAMPLES/${sample_id}"
    TRIM_DIR="${SAMPLE_DIR}/trim"
    MAP_DIR="${SAMPLE_DIR}/alignment"
    IN1="${sample_id}_1.fastq"
    IN2="${sample_id}_2.fastq"

    # Basic structure
    mkdir -p "$SAMPLE_DIR"

    # Enter the sample folder
    cd "$SAMPLE_DIR" || { echo "Error: Unable to enter in $SAMPLE_DIR"; exit 1; }

    # Checking the raw inputs
    if [[ ! -f "$IN1" || ! -f "$IN2" ]]; then
        echo "⚠️  Raw files not found for $sample_id (IN1=$IN1, IN2=$IN2)"
        continue
    fi
    echo "IN1: $IN1"
    echo "IN2: $IN2"
    echo "STRAND: [$STRAND]"

    # Performing trimming only if trimmed fastq files do not exist
    if [[ ! -f "$TRIM_DIR/${sample_id}_trimmed_1.fastq" || ! -f "$TRIM_DIR/${sample_id}_trimmed_2.fastq" ]]; then
        mkdir -p "$TRIM_DIR" "$MAP_DIR"

        # ------------ 1) TRIMMING READS ---------------------------
        case "$STRAND" in
          forward)
            bash "$TRIM_FW_SCRIPT" "$sample_id" "$IN1" "$IN2" "$TRIM_DIR" "$adapt_r1" "$adapt_r2" "$trimming_type" "$threads" "$trimming_quality_threshold" "$min_length_trim"
            ;;
          reverse)
            bash "$TRIM_RV_SCRIPT" "$sample_id" "$IN1" "$IN2" "$TRIM_DIR" "$adapt_r1" "$adapt_r2" "$trimming_type" "$threads" "$trimming_quality_threshold" "$min_length_trim"
            ;;
          *)
            echo "[WARN] $sample_id: unknown strandedness '$STRAND' (expected forward/reverse)" >&2
            continue
            ;;
        esac
    else
        echo "✓ Trimmed FASTQs already exist for $sample_id — skipping trimming."
    fi

done


# Quality control before alignment

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
        bash "$MAP_SCRIPT" "$sample_id" "${TRIM_DIR}/${sample_id}_trimmed_1.fastq" "${TRIM_DIR}/${sample_id}_trimmed_2.fastq" "$threads" "$MAP_DIR" "$GENOME_DIR" "$mismatch_align"
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
        bash "$MARKDUP_SCRIPT" "$sample_id" "$threads" "$MAP_DIR" "$remove_duplicates"
  fi
  done


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
   else
   echo 'Skipping quantification step'
   fi

  done



# ---------- 5) COUNT MATRIXES MERGE ----------------------

   if [ ! -d "$CWD/SAMPLES/DEA_results" ]; then

      mkdir $CWD/SAMPLES/DEA_results
   fi

   DEA_results="$CWD/SAMPLES/DEA_results"
      
   REP_GTF_PATH=$CWD/$repeat_gtf 

   if [ ! -f "$CWD/SAMPLES/count_data.txt" ]; then
     echo 'count_data.txt must be generated'

     Rscript "$MERGE_QUANT_SCRIPT" "$CWD" "$SAMPLES_DIR" "$REP_GTF_PATH" "$threads" "$DEA_results" "$sample_list"
     cd ..
   else
   echo 'count_data.txt already exists'
   fi



# ---------- 6) DIFFERENTIAL EXPRESSION ANALYSIS ----------

cond=$(awk '
BEGIN { FS = "\t" }  # Usa tabulador; cambia a FS="," si es CSV
NR==1 {
  for (i=1; i<=NF; i++) {
    if ($i == "condition") col=i
  }
  next
}
{
  if (col) print $col
}
' $CWD/$sample_list | sort | uniq | wc -l)

if [ "$cond" -eq 2 ]; then

Rscript "$DEA_SCRIPT" "$batch_effect" "$sample_list" "$method" "$CWD" "$REP_GTF_PATH" "$k_num" "$FDR" "$log2FC"

else
Rscript "$DEA_SCRIPT_multicond" "$batch_effect" "$sample_list" "$method" "$CWD" "$REP_GTF_PATH" "$k_num" "$FDR" "$log2FC"
fi



# ---------- 7) TRANSCRIPTOME ASSEMBLY ---------------------


cat "$CWD/$genome_gtf" "$CWD/$repeat_gtf" > "$CWD/combined_annotations.gtf"

combined_annotations="$CWD/combined_annotations.gtf"

if [ ! -d "$CWD/Transcriptome_assembly" ]; then

mkdir $CWD/Transcriptome_assembly

fi
 
 awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1, $2, $3}' "../$sample_list" \
  | while IFS=$'\t' read -r sample_id STRAND CONDITION; do
      [[ -z "$sample_id" ]] && continue

      SAMPLE_DIR="$CWD/SAMPLES/${sample_id}"

      if [[ ! -f "$SAMPLE_DIR/${sample_id}_transcriptome.gtf" ]];then

          bash "$ASSEMBLY_SCRIPT" "$combined_annotations" "$sample_id" "$threads" "$STRAND"
          
          cp $SAMPLE_DIR/${sample_id}_transcriptome.gtf $CWD/Transcriptome_assembly
      
      fi
done


if [ ! -d "$CWD/Transcriptome_assembly_novels" ]; then

mkdir $CWD/Transcriptome_assembly_novels

fi

bash "$ASSEMBLY_SCRIPT_2" "$threads" "$PREPDE_SCRIPT" "$CWD" "$genome_gtf" "$repeat_gtf" "$sample_list" "$STRAND" 





# ---------- 7) WGCNA COEXPRESSION ANALYSIS ---------------

DEA_DIR=$SAMPLES_DIR/DEA_results

if [ ! -d "$CWD/Coexpression_analysis" ]; then

mkdir $CWD/Coexpression_analysis

fi

COEXPRESSION_DIR=$CWD/Coexpression_analysis

Rscript "$WGCNA_SCRIPT" "$DEA_DIR" "$CWD" "$sample_list" "$COEXPRESSION_DIR"




# ---------- 8) PREDICTION MODEL ANALYSIS -----------------

if [ "$prediction_model" = "yes" ]; then
rows=$(( $(wc -l < "$CWD/$sample_list") - 1 ))
if [ "$rows" -gt 40 ]; then
Rscript "$PRED_MODEL" "$CWD" "$sample_list" "$threads"
fi
fi

}


#Running the function

run_pipeline_sample 
