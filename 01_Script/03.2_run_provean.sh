#!/bin/bash

# ==============================================================================
# 03.2_run_provean.sh - VERSION FINALE PROPRE
# ==============================================================================

CONDA_BASE=$(conda info --base 2>/dev/null)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate /data/projet1/conda/provean_env

VAR_DIR="/data/projet1/02_Pipeline_Output/03_Mutations"
REF_FASTA="/data/projet1/00_Data/ref_seq_proteique_fasta/orf_trans_all.fasta"
BLAST_DB="/data/projet1/00_Data/ref_seq_proteique_fasta/Ref_S288C_Prot"
OUT_DIR="/data/projet1/02_Pipeline_Output/04_Provean"
THREADS=4

mkdir -p "$OUT_DIR"
TMP_DIR="/tmp/provean_run_$$"
mkdir -p "$TMP_DIR"
trap "rm -rf $TMP_DIR" EXIT

echo "=================================================="
echo "DÉMARRAGE PROVEAN"
echo "=================================================="

for var_file in "$VAR_DIR"/*.var; do
    [ -e "$var_file" ] || continue
    
    filename=$(basename "$var_file")
    strain_name="${filename%.var}"
    out_file="$OUT_DIR/${strain_name}_provean.tsv"
    log_file="$OUT_DIR/${strain_name}.log"
    
    echo "-> Traitement de la souche : $strain_name"
    > "$out_file"
    > "$log_file"
    
    clean_var="$TMP_DIR/${strain_name}_clean.var"
    tr -d '\r' < "$var_file" > "$clean_var"
    
    prot_ids=$(cut -f1 "$clean_var" | sort -u)
    
    for prot_id in $prot_ids; do
        [ -z "$prot_id" ] && continue
        
        single_fasta="$TMP_DIR/${prot_id}.fasta"
        single_var="$TMP_DIR/${prot_id}.var"
        
        echo ">${prot_id}" > "$single_fasta"
        awk -v id="$prot_id" 'BEGIN{p=0} /^>/{if($0 ~ ">"id"($|[ \t])") {p=1} else {p=0}} {if(p && !/^>/) print}' "$REF_FASTA" | tr -d '\r' >> "$single_fasta"
        
        grep "^${prot_id}[[:space:]]" "$clean_var" | cut -f2 > "$single_var"
        
        if [ -s "$single_fasta" ] && [ -s "$single_var" ]; then
            provean -q "$single_fasta" -d "$BLAST_DB" -v "$single_var" --num_threads "$THREADS" >> "$out_file" 2>> "$log_file"
        fi
    done
done

echo "=================================================="
echo "PIPELINE TERMINÉ."
echo "=================================================="