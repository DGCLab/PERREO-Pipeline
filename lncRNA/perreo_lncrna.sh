#!/usr/bin/env bash
set -euo pipefail

#Defining current path
CWD="$(pwd)"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/scripts_PE/logging.sh"

# Scripts paths that already exist (ajusta nombres/paths)
TRIM_FW_SCRIPT="$CWD/scripts_PE/script_trimming_forward.sh"   
TRIM_RV_SCRIPT="$CWD/scripts_PE/script_trimming_reverse.sh"
TRIM_UNSTR_SCRIPT="$CWD/scripts_PE/script_trimming_unstranded.sh"
TRIM_EXTRA="$CWD/scripts_PE/trimGC.py"
MAP_SCRIPT="$CWD/scripts_PE/script_map_pe.sh"
MARKDUP_SCRIPT="$CWD/scripts_PE/script_markduplicates.sh"           
QUANT_SCRIPT="$CWD/scripts_PE/quant.R"
MERGE_QUANT_SCRIPT="$CWD/scripts_PE/merge_quant.R"
DEA_SCRIPT="$CWD/scripts_PE/dea.R"
DEA_SCRIPT_multicond="$CWD/scripts_PE/dea_multicond.R"
ASSEMBLY_SCRIPT="$CWD/scripts_PE/stringtie2.sh"   
WGCNA_SCRIPT="$CWD/scripts_PE/WGCNA.R"         
PRED_MODEL="$CWD/scripts_PE/prediction_model.R"
HYBRIDS_SCRIPT="$CWD/scripts_PE/hybrid_transcripts.sh"
HYBRIDS_R_SCRIPT="$CWD/scripts_PE/hybrid_transcripts_visualization_script.R"

# Parsing arguments
while [[ $# -gt 0 ]]; do
  case $1 in
      -sample_list) sample_list="$2"; shift 2 ;;
      -reference_genome) reference_genome="$2"; shift 2 ;;
      -genome_gtf) genome_gtf="$2"; shift 2 ;;
      -lnc_gtf) lnc_gtf="$2"; shift 2 ;;
      -threads) threads="$2"; shift 2 ;;
      -ram) threads="$2"; shift 2 ;;
      -adapt_r1) adapt_r1="$2"; shift 2 ;;
      -adapt_r2) adapt_r2="$2"; shift 2 ;;
      -trimming_quality_threshold) trimming_quality_threshold="$2"; shift 2 ;;
      -min_length_trim) min_length_trim="$2"; shift 2 ;;
      -max_length_trim) min_length_trim="$2"; shift 2 ;;
      -initial_trim_read1) initial_trim_read1="$2"; shift 2 ;;
      -initial_trim_read2) initial_trim_read2="$2"; shift 2 ;;
      -trimming) trimming="$2"; shift 2 ;;
      -mismatch_align) mismatch_align="$2"; shift 2 ;;
      -project_name) project_name="$2"; shift 2 ;;
      -remove_duplicates) remove_duplicates="$2"; shift 2 ;;
      -log2FC) log2FC="$2"; shift 2 ;;
      -FDR) FDR="$2"; shift 2 ;;
      -batch) batch="$2"; shift 2 ;;
      -method) method="$2"; shift 2 ;;
      -prediction_model) prediction_model="$2"; shift 2 ;;
      -positive_class) positive_class="$2"; shift 2 ;;
      -polya) polya="--poly-a"; shift 2 ;;
      *) echo "Unknown argument: $1"; shift ;;
  esac
done


# --- Assigning values by default ---
threads=${threads:-8}
ram=${ram:-32}
k_num=${k_num:-2}
FDR=${FDR:-0.05}
log2FC=${log2FC:-1.0}
mismatch_align=${mismatch_align:-0.05}
trimming_quality_threshold=${trimming_quality_threshold:-30}
min_length_trim=${min_length_trim:-16}
initial_trim_read1=${initial_trim_read1:-0}
initial_trim_read2=${initial_trim_read2:-0}
prediction_model=${prediction_model:-no}
batch=${batch:-no}
trimming=${trimming:-simple}
polya="${polya:-}"
positive_class="${positive_class:-}"
adapt_r1="${adapt_r1:-}"
adapt_r2="${adapt_r2:-}"

if [[ "$prediction_model" == "yes" && -z "$positive_class" ]]; then
  echo "ERROR: --positive_class is mandatory when --prediction_model yes" >&2
  exit 2
fi


# Function for running the analysis with all the samples

msg_info "Starting PERREO pipeline"

run_pipeline_sample() {
  msg_info " ----- Pipeline paramters summary -----"
  msg_info "samplesheet: ${GREEN}$sample_list${RESET}"
  msg_info "Reference: ${GREEN}$reference_genome${RESET}"
  msg_info "GTF: ${GREEN}$genome_gtf${RESET}"
  msg_info "Repeats: ${GREEN}$lnc_gtf${RESET}"
  msg_info "Threads: ${GREEN}${threads-8}${RESET}"
  msg_info "R1 adaptor: ${GREEN}${adapt_r1:-none}${RESET}"
  msg_info "R2 adaptor: ${GREEN}${adapt_r2:-none}${RESET}"
  msg_info "trimming type: ${GREEN}${trimming:-simple}${RESET}"
  msg_info "trimming_quality_threshold: ${GREEN}${trimming_quality_threshold:-30}${RESET}" 
  msg_info "minimum_length_trim: ${GREEN}${min_length_trim:-16}${RESET}"
  msg_info "initial_trim_read1: ${GREEN}${initial_trim_read1:-0}${RESET}"
  msg_info "initial_trim_read2: ${GREEN}${initial_trim_read2:-0}${RESET}"
  msg_info "mismatch_align type: ${GREEN}${mismatch_align:-0.05}${RESET}"
  msg_info "Project: ${GREEN}$project_name${RESET}"
  msg_info "Remove duplicates: ${GREEN}$remove_duplicates${RESET}"
  msg_info "Batch effect: ${GREEN}${batch:-no}${RESET}"
  msg_info "Method: ${GREEN}$method${RESET}"
  msg_info "log2FC: ${GREEN}${log2FC:-1}${RESET}"
  msg_info "FDR: ${GREEN}${FDR:-0.05}${RESET}"
  msg_info "Positive class: ${GREEN}$positive_class${RESET}"

###

if [ -d "genome_index" ]; then
  msg_ok "Genome index already exists in genome_index, creation omitted."
else
  msg_info "[STAR] Creating genome index..."
  STAR --runMode genomeGenerate --runThreadN $threads --genomeDir genome_index --genomeFastaFiles "$reference_genome" --sjdbGTFfile "$genome_gtf" --sjdbOverhang 100
fi

GENOME_DIR=$(realpath genome_index)

cd samples
SAMPLES_DIR=$CWD/samples

set -euo pipefail

awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1, $2, $3}' "../$sample_list" \
| while IFS=$'\t' read -r sample_id STRAND CONDITION; do
    [[ -z "$sample_id" ]] && continue

    cutadapt_threads=$threads

    STRAND=$(echo "$STRAND" | tr '[:upper:]' '[:lower:]' | xargs)
    CONDITION=$(echo "${CONDITION:-}" | xargs)

    msg_info "[CUTADAPT] Processing ${GREEN}$sample_id${RESET} with strandedness=${GREEN}$STRAND${RESET} and condition=${GREEN}$CONDITION${RESET}"

    # Path per sample
    SAMPLE_DIR="$CWD/samples/${sample_id}"
    TRIM_DIR="${SAMPLE_DIR}/trim"
    MAP_DIR="${SAMPLE_DIR}/alignment"
    IN1="${sample_id}_1.fastq"
    IN2="${sample_id}_2.fastq"

    # Basic structure
    mkdir -p "$SAMPLE_DIR"

    # Enter the sample folder
    cd "$SAMPLE_DIR" || { msg_error "[CUTADAPT] Unable to enter in $SAMPLE_DIR"; exit 1; }

    # Checking the raw inputs
    if [[ ! -f "$IN1" || ! -f "$IN2" ]]; then
        msg_error "[CUTADAPT] Raw files not found for $sample_id (IN1=$IN1, IN2=$IN2)"
        continue
    fi

    # Performing trimming only if trimmed fastq files do not exist
    if [[ ! -f "$TRIM_DIR/${sample_id}_trimmed_1.fastq" || ! -f "$TRIM_DIR/${sample_id}_trimmed_2.fastq" ]]; then
        mkdir -p "$TRIM_DIR" "$MAP_DIR"

        # ------------ 1) TRIMMING READS ---------------------------
        case "$STRAND" in
          forward)
            bash "$TRIM_FW_SCRIPT" "$sample_id" "$IN1" "$IN2" "$TRIM_DIR" "$adapt_r1" "$adapt_r2" "$trimming" "$threads" "$trimming_quality_threshold" "$min_length_trim" "$initial_trim_read1" "$initial_trim_read2" "$polya"
            ;;
          reverse)
            bash "$TRIM_RV_SCRIPT" "$sample_id" "$IN1" "$IN2" "$TRIM_DIR" "$adapt_r1" "$adapt_r2" "$trimming" "$threads" "$trimming_quality_threshold" "$min_length_trim" "$initial_trim_read1" "$initial_trim_read2" "$polya"
            ;;
            unstranded)
            bash "$TRIM_UNSTR_SCRIPT" "$sample_id" "$IN1" "$IN2" "$TRIM_DIR" "$adapt_r1" "$adapt_r2" "$threads" "$trimming_quality_threshold" "$min_length_trim" "$initial_trim_read1" "$initial_trim_read2" "$polya"
            ;;
          *)
            msg_warn "[CUTADAPT] $sample_id: unknown strandedness '$STRAND' (expected forward/reverse)" >&2
            continue
            ;;
        esac
    else
        msg_ok "[CUTADAPT] Trimmed FASTQs already exist for $sample_id — skipping trimming."
    fi

done


# Quality control before alignment

   cd $SAMPLES_DIR


# ---------- 2) MAPPING AGAINST REFERENCE GENOME ----------

msg_info "[STAR] Starting alignment against reference genome"

  awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1, $2, $3}' "../$sample_list" \
  | while IFS=$'\t' read -r sample_id STRAND CONDITION; do
     [[ -z "$sample_id" ]] && continue

     STRAND=$(msg_info "[STAR] $STRAND" | tr '[:upper:]' '[:lower:]' | xargs)
     SAMPLE_DIR="$CWD/samples/${sample_id}"
     TRIM_DIR="${SAMPLE_DIR}/trim"
     MAP_DIR="${SAMPLE_DIR}/alignment"

  if [[ -f "$MAP_DIR/${sample_id}_Aligned.sortedByCoord.out.bam" ]];then
        msg_info '[STAR] Skipping mapping...'
  else
        bash "$MAP_SCRIPT" "$sample_id" "${TRIM_DIR}/${sample_id}_trimmed_1.fastq" "${TRIM_DIR}/${sample_id}_trimmed_2.fastq" "$threads" "$MAP_DIR" "$GENOME_DIR" "$mismatch_align"
  fi
  done

  msg_ok "[STAR] Alignment for all samples completed"
  
# -----------3) MARKDUPLICATES ----------------------------

  
 awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1, $2, $3}' "../$sample_list" \
  | while IFS=$'\t' read -r sample_id STRAND CONDITION; do
     [[ -z "$sample_id" ]] && continue
  
  SAMPLE_DIR="$CWD/samples/${sample_id}"
  MAP_DIR="${SAMPLE_DIR}/alignment"

  if [[ -f "$MAP_DIR/${sample_id}_marked_duplicates_STAR.bam" ]];then
        msg_info '[MARKDUPLICATES] Skipping Markduplicates...'
  else
        REMOVE_DUPLICATES=false
        bash "$MARKDUP_SCRIPT" "$sample_id" "$threads" "$MAP_DIR" "$remove_duplicates" "$ram"
  fi
  done


# ---------- 4) QUANTIFICATION ----------------------------

  awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1, $2, $3}' "../$sample_list" \
  | while IFS=$'\t' read -r sample_id STRAND CONDITION; do
      [[ -z "$sample_id" ]] && continue

      SAMPLE_DIR="$CWD/samples/${sample_id}"
      MAP_DIR="${SAMPLE_DIR}/alignment"

   if [ ! -d "${SAMPLE_DIR}/Quantification" ]; then 
      mkdir ${SAMPLE_DIR}/Quantification
   fi
   
      msg_info "[FEATURECOUNTS] Starting ${sample_id} quantification..."
      
      QUANT_DIR="${SAMPLE_DIR}/Quantification"

      STRAND=$(echo "$STRAND" | tr '[:upper:]' '[:lower:]' | xargs)

      REP_GTF_PATH=$CWD/$lnc_gtf 

   if [ ! -f "$QUANT_DIR/${sample_id}_quant.txt" ]; then

   #Running Rsubread for quantification
      Rscript "$QUANT_SCRIPT" "$MAP_DIR" "$sample_id" "$REP_GTF_PATH" "$threads" "$STRAND" "$SAMPLE_DIR" "$QUANT_DIR"
   else
   msg_warn "[FEATURECOUNTS] ${sample_id}_quant.txt alredy exists, skipping."
   fi

  done



# ---------- 5) COUNT MATRIXES MERGE ----------------------

   if [ ! -d "$CWD/Results" ]; then
      mkdir $CWD/Results
   fi
   
   if [ ! -d "$CWD/Results/DEA_results" ]; then
      mkdir $CWD/Results/DEA_results
   fi

   DEA_results="$CWD/Results/DEA_results"
      
   REP_GTF_PATH=$CWD/$lnc_gtf 

   if [ ! -f "$CWD/Results/count_data.txt" ]; then
     msg_warn '[FEATURECOUNTS] count_data.txt must be generated'
     msg_info '[FEATURECOUNTS] generating...'
     Rscript "$MERGE_QUANT_SCRIPT" "$CWD" "$SAMPLES_DIR" "$REP_GTF_PATH" "$threads" "$sample_list"
     cd ..
   else
   msg_ok '[FEATURECOUNTS] count_data.txt already exists'
   fi

# ---------- 6) TRANSCRIPTOME ASSEMBLY ---------------------

if [ ! -d "$CWD/Results/transcriptome_assembly" ]; then
      mkdir $CWD/Results/transcriptome_assembly
   fi
   
cat "$CWD/$genome_gtf" "$CWD/$lnc_gtf" > "$CWD/Results/transcriptome_assembly/combined_annotations.gtf"

combined_annotations="$CWD/Results/transcriptome_assembly/combined_annotations.gtf"
 
 awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1, $2, $3}' "$CWD/$sample_list" \
  | while IFS=$'\t' read -r sample_id STRAND CONDITION; do
      [[ -z "$sample_id" ]] && continue

      SAMPLE_DIR="$CWD/samples/${sample_id}"

      if [[ ! -f "$CWD/Results/transcriptome_assembly/${sample_id}_transcriptome.gtf" ]];then
          
          msg_info "[STRINGTIE2] Starting transcriptome assembly..."
          bash "$ASSEMBLY_SCRIPT" "$combined_annotations" "$sample_id" "$threads" "$STRAND" "$CWD"
          
          mv $SAMPLE_DIR/${sample_id}_transcriptome.gtf $CWD/Results/transcriptome_assembly
      
      fi
done

msg_ok "[STRINGTIE2] All .gtf generated, transcriptome assembly was generated successfully."


#Then, we call this script to quantify reads that map uniquely to genes, uniquely to repeats and reads that map to both genes and repetitive regions

if [[ ! -f "$CWD/Results/transcriptome_assembly/hybrid_transcripts_summary.tsv" ]]; then

msg_info "Calculating hybrid transcripts..."
bash "$HYBRIDS_SCRIPT" "$CWD/Results/transcriptome_assembly" "$CWD" "$genome_gtf" "$lnc_gtf" > $CWD/Results/transcriptome_assembly/hybrid_transcripts_summary.tsv

fi

Rscript "$HYBRIDS_R_SCRIPT" "$CWD"
msg_ok "Done!"


# ---------- 7) DIFFERENTIAL EXPRESSION ANALYSIS ----------

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

Rscript "$DEA_SCRIPT" "$batch" "$sample_list" "$method" "$CWD" "$REP_GTF_PATH" "$k_num" "$FDR" "$log2FC"

else
Rscript "$DEA_SCRIPT_multicond" "$batch" "$sample_list" "$method" "$CWD" "$REP_GTF_PATH" "$k_num" "$FDR" "$log2FC"
fi

msg_ok "DEA Analysis completed"



# ---------- 8) WGCNA COEXPRESSION ANALYSIS ---------------

DEA_DIR=$CWD/Results/DEA_results

msg_info "Starting WGCNA coexpression analysis..."

if [ ! -d "$CWD/Results/Coexpression_analysis" ]; then

mkdir -p $CWD/Results/Coexpression_analysis

fi

COEXPRESSION_DIR=$CWD/Results/Coexpression_analysis

Rscript "$WGCNA_SCRIPT" "$DEA_DIR" "$CWD" "$sample_list" "$COEXPRESSION_DIR"

msg_ok "WGCNA coexpression analysis completed"


# ---------- 9) PREDICTION MODEL ANALYSIS -----------------

if [ "$prediction_model" = "yes" ]; then
rows=$(( $(wc -l < "$CWD/$sample_list") - 1 ))
if [ "$rows" -gt 40 ]; then
msg_info "Starting prediction model analysis..."
Rscript "$PRED_MODEL" "$CWD" "$sample_list" "$threads" "$positive_class"

msg_ok "Prediction model generated"
fi
fi



msg_ok "PERREO successfully completed"


}


#Running the function

run_pipeline_sample 
