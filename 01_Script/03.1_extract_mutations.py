#!/usr/bin/env python3

#input = BLASTp output tsv file
#output = .var file per strains 


import os 
import argparse 
import glob
from Bio import SeqIO
from Bio import Align



def extract_mutations(ref_aln, query_aln):
    mutations = []
    ref_pos = 0  
    
    for r, q in zip(ref_aln, query_aln):
        if r != '-':
            ref_pos += 1
        
        if r == '*' or q == '*':
            continue
            
        if r != q:  
            if r == '-': # insertion
                mutations.append(f"{ref_pos}_{ref_pos+1}ins{q}")
            elif q == '-': # deletion
                mutations.append(f"{r}{ref_pos}del")
            else: # substitution
                mutations.append(f"{r}{ref_pos}{q}")
                
    return mutations

def process_souche(blastp_tsv,souche_fasta, ref_dico, out_var,aligner):
    #traite une souche, lit tsv, aligne seq et crée fichier .var

    try:
        strain_dico = SeqIO.to_dict(SeqIO.parse(souche_fasta, "fasta"))
    except FileNotFoundError:
        print(f"Fichier FASTA introuvable pour la souche : {souche_fasta}")
        return
    mutations_count = 0

    with open(blastp_tsv,'r') as f_tsv, open (out_var,'w') as f_out:
        for line in f_tsv:
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue  # ignire les lignes mal formées
            qseqid = parts[0]
            sseqid = parts[1]

            if qseqid in strain_dico and sseqid in ref_dico:
                ref_seq = str(ref_dico[sseqid].seq)
                query_seq = str(strain_dico[qseqid].seq)

                if ref_seq == query_seq:
                    continue #aucune mutation à déclarer, veuillez circuler (je rigole)

                alignements = aligner.align(ref_seq, query_seq)
                best_alignment = alignements[0]

                ref_alignement = best_alignment[0]
                query_alignement = best_alignment[1]

                mutations = extract_mutations(ref_alignement, query_alignement)

                for mutants in mutations: 
                    f_out.write(f"{sseqid}\t{mutants}\n")
                    mutations_count += 1

    print(f"  -> {mutations_count} mutations trouvées et enregistrées.")




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extrait les mutations à partir des fichiers TSV de BLASTp et génère des fichiers .var pour chaque souche.")

    parser.add_argument("-t", "--tsv_dir", default="/data/projet1/02_Pipeline_Output/02_Blast", help="Dossier contenant les fichiers TSV de BLASTp (défaut: %(default)s)")
    parser.add_argument("-q", "--query_dir", default="/data/projet1/02_Pipeline_Output/traduction_test_batch", help="Dossier des FASTA protéiques des souches")
    parser.add_argument("-r", "--ref_fasta", default="/data/projet1/00_Data/reference_S288C/orf_trans_all.fasta", help="Fichier FASTA de référence S288C (défaut: %(default)s)")
    parser.add_argument("-o", "--out_dir", default="/data/projet1/02_Pipeline_Output/03_Mutations", help="Dossier de sortie pour les fichiers .var (défaut: %(default)s)")

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("Chargement de la séquence de référence...")
    ref_dico = SeqIO.to_dict(SeqIO.parse(args.ref_fasta, "fasta"))

    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -10     #surement besoin de modifier ces paramètres 
    aligner.extend_gap_score = -0.5  # en fonction de ce qui se passe
    
    tsv_files = glob.glob(os.path.join(args.tsv_dir, "*.tsv"))

    for tsv_file in tsv_files:
        filename = os.path.basename(tsv_file)
        strain_name = filename.replace("_blast.tsv", "")
        strain_fasta = os.path.join(args.query_dir, f"{strain_name}.fasta")
        out_var = os.path.join(args.out_dir, f"{strain_name}.var")

        print(f"Traitement de la souche : {strain_name}")
        process_souche(tsv_file, strain_fasta, ref_dico, out_var, aligner)

    print("Etape 3 : Extraction des mutations terminée.")
