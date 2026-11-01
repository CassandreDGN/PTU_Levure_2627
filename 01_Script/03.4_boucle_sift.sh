#!/usr/bin/env bash

MUT_DIR="/data/projet1/02_Pipeline_Output/3_SIFT_Mutations"
OUT_DIR="/data/projet1/02_Pipeline_Output/4_SIFT_SWISSPROT_Results"
REF_FASTA="/data/projet1/00_Data/ref_seq_proteique_fasta/orf_trans_all.strict_headers.fasta"
DB_SIFT="/data/projet1/00_Data/swiss_prot/uniprot_sprot.fasta"

mkdir -p "$OUT_DIR"

export MUT_DIR OUT_DIR REF_FASTA DB_SIFT

process_strain() {
    strain_dir="$1"
    strain_name=$(basename "$strain_dir")

    # Correction du slash
    if [[ ! -d "${strain_dir}/subst" ]]; then
        return
    fi

    strain_out="$OUT_DIR/${strain_name}"

    if [[ -d "$strain_out" ]] && [[ $(ls -A "$strain_out" 2>/dev/null) ]]; then
        echo " [SKIP] ${strain_name} -> Déjà traité."
        return
    fi

    mkdir -p "$strain_out"
    echo "=== Processing SIFT4G for: ${strain_name} ==="
    
    # Correction du slash et sauvegarde des erreurs
    sift4g -q "$REF_FASTA" --subst "${strain_dir}/subst/" -d "$DB_SIFT" --out "$strain_out/" > /dev/null 2> "${OUT_DIR}/${strain_name}_sift_errors.log"
}

export -f process_strain

# -P 2 impose STRICTEMENT 2 cœurs en parallèle
find "$MUT_DIR" -mindepth 1 -maxdepth 1 -type d | xargs -P 1 -I {} bash -c 'process_strain "$@"' _ {}

echo "✅ SIFT4G avec Swiss-Prot terminé pour toutes les souches."