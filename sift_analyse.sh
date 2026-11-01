#!/bin/bash
# Script d'analyse SIFT4G pour le projet PTU
# Auteur : Beray SÖNMEZ
# Date : 28 octobre 2026

SIFT4G="/data/projet1/06_Conda/sift_env/bin/sift4g"
VAR_DIR="/data/projet1/02_Pipeline_Output/4_Mutations_par_genes"
OUT_DIR="/data/projet1/02_Pipeline_Output/4_SIFT_Results"

mkdir -p $OUT_DIR

if [ -z "$1" ]; then
    echo "Usage: $0 <souche>"
    exit 1
fi

SOUCHE=$1
GENE_DIR="$VAR_DIR/$SOUCHE.HP1.nuclear_genome"

for VAR_FILE in "$GENE_DIR"/*.var; do
    BASE=$(basename "$VAR_FILE" .var)
    echo "Traitement de $BASE..."
    $SIFT4G -q "$GENE_DIR/$BASE.fasta" -d /data/projet1/00_Data/db_sift/uniprot_sprot.fasta --subst "$GENE_DIR/$BASE.subst"
done

echo "Analyse terminée."
