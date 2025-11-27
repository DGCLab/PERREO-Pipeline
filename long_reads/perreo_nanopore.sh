#!/usr/bin/env bash
set -euo pipefail

#Defining current path
CWD="$(pwd)"
#GLOBAL_DIR="$CWD"

# Scripts paths that already exist (ajusta nombres/paths)
DEA_SCRIPT="$CWD/scripts_long_reads/dea.R"
#BAMBU_SCRIPT="$CWD/scripts_long_reads/script_bambu.R"
DEA_SCRIPT_multicond="$CWD/scripts_long_reads/dea_multicond.R"      
PRED_MODEL="$CWD/scripts_long_reads/prediction_model.R"     
WGCNA_SCRIPT="$CWD/scripts_long_reads/WGCNA.R"     
ASSEMBLY_SCRIPT="$CWD/scripts_long_reads/stringtie2.sh"   
QUANT_SCRIPT="$CWD/scripts_long_reads/quant.R"       
MERGE_QUANT_SCRIPT="$CWD/scripts_long_reads/merge_quant.R" 
      
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
      *) echo "Opción desconocida: $1"; shift ;;
  esac
done

# --- Asignar valores por defecto ---
threads=${threads:-8}
k_num=${k_num:-2}
FDR=${FDR:-0.05}
log2FC=${log2FC:-1.0}
prediction_model=${prediction_model:-yes}

# Function for running the analysis with all the samples

run_pipeline_sample() {
  echo " ----- Pipeline paramters summary -----"
  echo "samplesheet: $sample_list"
  echo "Reference: $reference_genome"
  echo "Repeats: $repeat_gtf"
  echo "Threads: $threads"
  echo "Proyect: $project_name"
  echo "Batch effect: $batch_effect"
  echo "Method: $method"
  echo "log2FC: $log2FC"
  echo "FDR: $FDR"

if [ -f "genome_index.mmi" ]; then
  echo "✅ Genome index already exists in genome_index, creation omitted."
else
  echo "Creating genome index..."
  minimap2 -d genome_index.mmi $reference_genome
fi


cd SAMPLES
SAMPLES_DIR=$CWD/SAMPLES


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

    SAMPLE_DIR="$CWD/SAMPLES/${sample_id}"

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

      SAMPLE_DIR="$CWD/SAMPLES/${sample_id}"

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


# ---------- 5) DIFFERENTIAL EXPRESSION ANALYSIS ----------


   if [ ! -d "$CWD/SAMPLES/DEA_results" ]; then

      mkdir $CWD/SAMPLES/DEA_results
   fi

   DEA_results="$CWD/SAMPLES/DEA_results"

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


# ---------- 6) TRANSCRIPTOME ASSEMBLY ------------------------


if [ ! -d "$CWD/Transcriptome_assembly" ]; then

mkdir $CWD/Transcriptome_assembly

fi

 awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1, $2, $3}' "$CWD/$sample_list" \
  | while IFS=$'\t' read -r sample_id STRAND CONDITION; do
      [[ -z "$sample_id" ]] && continue

      SAMPLE_DIR="$CWD/SAMPLES/${sample_id}"

      if [[ ! -f "$SAMPLE_DIR/${sample_id}_transcriptome.gtf" ]];then

          bash "$ASSEMBLY_SCRIPT" "$CWD" "$repeat_gtf" "$sample_id" "$threads" "$STRAND"
          
          cp $SAMPLE_DIR/${sample_id}_transcriptome.gtf $CWD/Transcriptome_assembly
      fi
    done







# ----------- 7) WGCNA ----------------------------------------

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
