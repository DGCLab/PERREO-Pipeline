###############################################################################
####                                                                       ####
####                      DIFFERENTIAL EXPRESSION ANALYSIS                 ####
####                                                                       ####
####             Authors: Francisco Rodríguez Martín                       ####                      
####                      Mario Masero León                                ####
####                                                                       ####
###############################################################################


#!/usr/bin/env bash
set -euo pipefail

#Defining current path
CWD="$(pwd)"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/scripts_long_reads/logging.sh"

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
      -prediction_model) prediction_model="$2"; shift 2 ;;
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
batch_effect="${batch_effect:-no}"

if [[ "$prediction_model" == "yes" && -z "$positive_class" ]]; then
  echo "ERROR: --positive_class is mandatory when --prediction_model yes" >&2
  exit 2
fi

# Function for running the analysis with all the samples

perreo_banner "LR" "Long-read sequencing technologies"

run_pipeline_sample() {

  perreo_summary "Pipeline parameters" \
    "Mode=LR" \
    "Samplesheet=${sample_list}" \
    "Reference genome=${reference_genome}" \
    "Repeat GTF=${repeat_gtf}" \
    "Threads=${threads}" \
    "Project=${project_name:-unnamed}" \
    "Batch effect=${batch_effect:-no}" \
    "DEA method=${method:-auto}" \
    "log2FC=${log2FC}" \
    "FDR=${FDR}" \
    "Prediction model=${prediction_model}" \
    "Positive class=${positive_class:-none}"

if [ -f "genome_index.mmi" ]; then
  msg_ok "Genome index already exists in genome_index.mmi, creation omitted."
else
  msg_info "[MINIMAP2] Creating genome index..."
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


perreo_stage 1 8 "Quality control"
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


perreo_stage 2 8 "Mapping against reference genome"
# ---------- 2) MAPPING AGAINST REFERENCE GENOME ----------

  awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1, $2, $3}' "../$sample_list" \
  | while IFS=$'\t' read -r sample_id STRAND CONDITION; do
      [[ -z "$sample_id" ]] && continue

    SAMPLE_DIR="$CWD/samples/${sample_id}"

  if [[ -f "$SAMPLE_DIR/${sample_id}.bam" ]];then
        msg_skip "[MINIMAP2] ${sample_id}.bam already exists, skipping mapping."
  else
        msg_info "[MINIMAP2] Performing alignment for ${sample_id}"
        minimap2 -t 14 -ax splice -uf -k14 -p 0.8 -N 100 "$CWD/genome_index.mmi"  "$SAMPLE_DIR/${sample_id}.fastq" > "$SAMPLE_DIR/${sample_id}.sam"

        samtools sort "$SAMPLE_DIR/${sample_id}.sam" -o "$SAMPLE_DIR/${sample_id}.bam"

        samtools index "$SAMPLE_DIR/${sample_id}.bam"

        NanoPlot --bam "$SAMPLE_DIR/${sample_id}.bam" -t $threads -o "$SAMPLE_DIR/qc_nanoplot_bam" 
  fi
  done



perreo_stage 3 8 "FeatureCounts quantification"
# ---------- 3) FEATURECOUNTS QUANTIFICATION ------------------

  awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1, $2, $3}' "../$sample_list" \
  | while IFS=$'\t' read -r sample_id STRAND CONDITION; do
      [[ -z "$sample_id" ]] && continue

      SAMPLE_DIR="$CWD/samples/${sample_id}"

      STRAND=$(echo "$STRAND" | tr '[:upper:]' '[:lower:]' | xargs)

      REP_GTF_PATH=$CWD/$repeat_gtf 

   if [ ! -f "$SAMPLE_DIR/${sample_id}_quant.txt" ]; then
   msg_info "[FEATURECOUNTS] Running quantification for ${sample_id}"

   #Running Rsubread for quantification
      Rscript "$QUANT_SCRIPT" "$sample_id" "$REP_GTF_PATH" "$threads" "$STRAND" "$SAMPLE_DIR" 
   else
   msg_skip "[FEATURECOUNTS] ${sample_id}_quant.txt already exists, skipping."
   fi

  done



perreo_stage 4 8 "Merge count matrices"
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


perreo_stage 5 8 "Transcriptome assembly"
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


perreo_stage 6 8 "Differential expression analysis"
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



perreo_stage 7 8 "WGCNA coexpression analysis"
# ----------- 7) WGCNA ----------------------------------------

DEA_DIR=$CWD/Results/DEA_results

msg_info "Starting WGCNA coexpression analysis..."

   if [ ! -d "$CWD/Results/Coexpression_analysis" ]; then
     mkdir $CWD/Results/Coexpression_analysis
   fi

COEXPRESSION_DIR=$CWD/Results/Coexpression_analysis

Rscript "$WGCNA_SCRIPT" "$DEA_DIR" "$CWD" "$sample_list" "$COEXPRESSION_DIR"
msg_ok "WGCNA coexpression analysis completed"



perreo_stage 8 8 "Prediction model analysis"
# ---------- 8) PREDICTION MODEL ANALYSIS -----------------


if [ "$prediction_model" = "yes" ]; then
rows=$(( $(wc -l < "$CWD/$sample_list") - 1 ))
if [ "$rows" -gt 40 ]; then
Rscript "$PRED_MODEL" "$CWD" "$sample_list" "$threads"

msg_ok "Prediction model generated"

fi
fi

perreo_finished "PERREO successfully completed"

}


#Running the function

run_pipeline_sample 
