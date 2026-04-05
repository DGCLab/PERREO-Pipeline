#!/usr/bin/env bash
set -euo pipefail

#Defining current path
CWD="$(pwd)"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/scripts_PE/logging.sh"

# Create .log
exec > >(tee -a "${SCRIPT_DIR}/pipeline.log") \
     2> >(tee -a "${SCRIPT_DIR}/pipeline.err" >&2)

# Scripts paths that already exist (ajusta nombres/paths)
DEA_SCRIPT="$CWD/scripts_long_reads/dea.R"
DEA_SCRIPT_multicond="$CWD/scripts_long_reads/dea_multicond.R"      
PRED_MODEL="$CWD/scripts_long_reads/prediction_model.R"     
WGCNA_SCRIPT="$CWD/scripts_long_reads/WGCNA.R"     
ASSEMBLY_SCRIPT="$CWD/scripts_long_reads/stringtie2.sh"   
QUANT_SCRIPT="$CWD/scripts_long_reads/quant.R"       
MERGE_QUANT_SCRIPT="$CWD/scripts_long_reads/merge_quant.R"
HYBRIDS_SCRIPT="$CWD/scripts_long_reads/hybrid_transcripts.sh"
HYBRIDS_R_SCRIPT="$CWD/scripts_long_reads/hybrid_transcripts_visualization_script.R"
      
  # Parsing arguments
while [[ $# -gt 0 ]]; do
  case $1 in
      -sample_list) sample_list="$2"; shift 2 ;;
      -reference_genome) reference_genome="$2"; shift 2 ;;
      -repeat_gtf) repeat_gtf="$2"; shift 2 ;;
      -threads) threads="$2"; shift 2 ;;
      -project_name) project_name="$2"; shift 2 ;;
      -log2FC) log2FC="$2"; shift 2 ;;
      -FDR) FDR="$2"; shift 2 ;;
      -batch_effect) batch_effect="$2"; shift 2 ;;
      -method) method="$2"; shift 2 ;;
      -prediction_model) method="$2"; shift 2 ;;
      -positive_class) positive_class="$2"; shift 2 ;;
      *) echo "Opción desconocida: $1"; shift ;;
  esac
done

# --- Asignar valores por defecto ---
threads="${threads:-8}"
k_num="${k_num:-2}"
FDR="${FDR:-0.05}"
log2FC="${log2FC:-1.0}"
prediction_model="${prediction_model:-no}"
positive_class="${positive_class:-}"

if [[ "$prediction_model" == "yes" && -z "$positive_class" ]]; then
  echo "ERROR: --positive_class is mandatory when --prediction_model yes" >&2
  exit 2
fi

# Function for running the analysis with all the samples

msg_info "Starting PERREO pipeline"

msg_info $' 
____   _____  ____   ____   _____   ___  
|  _ \\ | ____||  _ \\ |  _ \\ | ____| / _ \\ 
| |_) ||  _|  | |_) || |_) ||  _|  | | | |
|  __/ | |___ |  _ < |  _ < | |___ | |_| |
|_|    |_____||_| \\_\\|_| \\_\\|_____| \\___/ 
'

msg_info "🔥 Starting PERREO pipeline for data derived from long-reads sequencing technologies"

msg_info "☕ Please ensure coffee levels are above threshold."


run_pipeline_sample() {

  msg_info " ----- Pipeline paramters summary -----"
  msg_info "samplesheet: ${GREEN}$sample_list${RESET}"
  msg_info "Reference: ${GREEN}$reference_genome${RESET}"
  msg_info "Repeats: ${GREEN}$repeat_gtf${RESET}"
  msg_info "Threads: ${GREEN}${threads-8}${RESET}"
  msg_info "Project: ${GREEN}$project_name${RESET}"
  msg_info "Batch effect: ${GREEN}${batch:-no}${RESET}"
  msg_info "Method: ${GREEN}$method${RESET}"
  msg_info "log2FC: ${GREEN}${log2FC:-1}${RESET}"
  msg_info "FDR: ${GREEN}${FDR:-0.05}${RESET}"
  msg_info "Positive class: ${GREEN}$positive_class${RESET}"

if [ -f "genome_index.mmi" ]; then
  echo "✅ Genome index already exists in genome_index, creation omitted."
else
  echo "Creating genome index..."
  minimap2 -d genome_index.mmi $reference_genome
fi


cd samples
SAMPLES_DIR=$CWD/samples


# ---------- 0) DEFINING VARIABLES ------------------------
awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1, $2, $3}' "../$sample_list" \
| while IFS=$'\t' read -r sample_id STRAND CONDITION; do
    [[ -z "$sample_id" ]] && continue

    STRAND=$(echo "$STRAND" | tr '[:upper:]' '[:lower:]' | xargs)
    CONDITION=$(echo "${CONDITION:-}" | xargs)
done


# ---------- 1) QUALITY CONTROL ---------------------------

  #awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1, $2, $3}' "../$sample_list" \
  #| while IFS=$'\t' read -r sample_id STRAND CONDITION; do
    #  [[ -z "$sample_id" ]] && continue

  # if [ ! -f "$SAMPLES_DIR/${sample_id}/QC_nanoplot" ];then

   #echo 'Running QC on ${sample_id}'

   #NanoPlot --fastq "$SAMPLES_DIR/${sample_id}/${sample_id}.fastq" -t $threads -o "$SAMPLES_DIR/${sample_id}/QC_nanoplot"
   
   #else
   #echo 'QC analyses already done'
   #fi

  #done


# ---------- 2) MAPPING AGAINST REFERENCE GENOME ----------

  awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1, $2, $3}' "../$sample_list" \
  | while IFS=$'\t' read -r sample_id STRAND CONDITION; do
      [[ -z "$sample_id" ]] && continue

    SAMPLE_DIR="$CWD/samples/${sample_id}"

  if [[ -f "$SAMPLE_DIR/${sample_id}.bam" ]];then
        echo 'Skipping mapping...'
  else
        echo "Performing alignment for $sample_id"
        minimap2 -t 14 -ax splice -uf -k14 -p 0.8 -N 100 "$CWD/genome_index.mmi"  "$SAMPLE_DIR/${sample_id}.fastq" > "$SAMPLE_DIR/${sample_id}.sam"

        samtools sort "$SAMPLE_DIR/${sample_id}.sam" -o "$SAMPLE_DIR/${sample_id}.bam"

        samtools index "$SAMPLE_DIR/${sample_id}.bam"

        NanoPlot --bam "$SAMPLE_DIR/${sample_id}.bam" -t $threads -o "$SAMPLE_DIR/qc_nanoplot_bam" 
  fi
  done



# ---------- 3) FEATURECOUNTS QUANTIFICATION ------------------

  awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1, $2, $3}' "../$sample_list" \
  | while IFS=$'\t' read -r sample_id STRAND CONDITION; do
      [[ -z "$sample_id" ]] && continue

      SAMPLE_DIR="$CWD/samples/${sample_id}"

      STRAND=$(echo "$STRAND" | tr '[:upper:]' '[:lower:]' | xargs)

      REP_GTF_PATH=$CWD/$repeat_gtf 

   if [ ! -f "$SAMPLE_DIR/${sample_id}_quant.txt" ]; then
   echo "Running quantification step"

   #Running Rsubread for quantification
      Rscript "$QUANT_SCRIPT" "$sample_id" "$REP_GTF_PATH" "$threads" "$STRAND" "$SAMPLE_DIR" 
   else
   echo 'Skipping quantification step'
   fi

  done



# -------- 4) QUANTIFICATION MERGE ---------------------------
#Finally, we call the last script to merge all the count matrixes

   if [ ! -d "$CWD/Results/DEA_results" ]; then

      mkdir $CWD/Results/DEA_results
   fi

   DEA_results="$CWD/Results/DEA_results"
      
   REP_GTF_PATH=$CWD/$repeat_gtf 

   if [ ! -f "$CWD/Results/count_data.txt" ]; then
     msg_warn '[FEATURECOUNTS] count_data.txt must be generated'
     msg_info '[FEATURECOUNTS] generating...'
     Rscript "$MERGE_QUANT_SCRIPT" "$CWD" "$SAMPLES_DIR" "$REP_GTF_PATH" "$threads" "$DEA_results" "$sample_list"
     cd ..
   else
   msg_ok '[FEATURECOUNTS] count_data.txt already exists'
   fi


# ---------- 5) TRANSCRIPTOME ASSEMBLY ------------------------

if [ ! -d "$CWD/Results/transcriptome_assembly" ]; then
      mkdir $CWD/Results/transcriptome_assembly
   fi
   
cat "$CWD/$genome_gtf" "$CWD/$repeat_gtf" > "$CWD/Results/transcriptome_assembly/combined_annotations.gtf"

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


# ---------- 6) DIFFERENTIAL EXPRESSION ANALYSIS ----------


   if [ ! -d "$CWD/Results/DEA_results" ]; then

      mkdir $CWD/Results/DEA_results
   fi

   DEA_results="$CWD/Results/DEA_results"

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



# ----------- 7) WGCNA ----------------------------------------

DEA_DIR=$CWD/Results/DEA_results

msg_info "Starting WGCNA coexpression analysis..."

   if [ ! -d "$CWD/Results/Coexpression_analysis" ]; then
     mkdir $CWD/Results/Coexpression_analysis
   fi

COEXPRESSION_DIR=$CWD/Results/Coexpression_analysis

Rscript "$WGCNA_SCRIPT" "$DEA_DIR" "$CWD" "$sample_list" "$COEXPRESSION_DIR"
msg_ok "WGCNA coexpression analysis completed"



# ---------- 8) PREDICTION MODEL ANALYSIS -----------------


if [ "$prediction_model" = "yes" ]; then
rows=$(( $(wc -l < "$CWD/$sample_list") - 1 ))
if [ "$rows" -gt 40 ]; then
Rscript "$PRED_MODEL" "$CWD" "$sample_list" "$threads"

msg_ok "Prediction model generated"

fi
fi

msg_ok "PERREO successfully completed"

}


#Running the function

run_pipeline_sample 
