#!/usr/bin/env python3
import os
import glob
import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

def run_blastp_single(args_tuple):
    strain_fasta, ref_db, output_tsv, blastp_bin, threads = args_tuple
    filename = os.path.basename(strain_fasta)
    
    cmd = [
        blastp_bin,
        "-query", strain_fasta,
        "-db", ref_db,
        "-out", output_tsv,
        "-outfmt", "6 qseqid sseqid pident length evalue bitscore qseq sseq",
        "-max_target_seqs", "1",
        "-evalue", "1e-5",
        "-num_threads", str(threads)
    ]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return filename, True
    except subprocess.CalledProcessError:
        return filename, False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Exécute le BLASTp en parallèle.")
    
    parser.add_argument("-i", "--data_dir", default="/data/projet1/02_Pipeline_Output/1_Translated_Proteins")
    parser.add_argument("-db", "--ref_db", default="/data/projet1/00_Data/ref_seq_proteique_fasta/Ref_S288C_Prot")
    parser.add_argument("-o", "--out_dir", default="/data/projet1/02_Pipeline_Output/2_Blast")
    parser.add_argument("-b", "--blastp_bin", default="/home/degen/.conda/envs/blastp_env/bin/blastp")
    parser.add_argument("-j", "--jobs", type=int, default=2)
    parser.add_argument("-t", "--threads", type=int, default=1)

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # Recherche directe des fichiers .prot.fasta
    fasta_files = sorted(glob.glob(os.path.join(args.data_dir, "*.prot.fasta")))

    if not fasta_files:
        print(f"Erreur : Aucun fichier .prot.fasta trouvé dans {args.data_dir}")
        exit(1)

    print(f"Lancement de BLASTp sur {len(fasta_files)} fichiers .prot.fasta...")
    print(f"Configuration : {args.jobs} jobs en parallèle, {args.threads} thread par job (2 cœurs au total).\n")

    tasks = []
    for fasta_file in fasta_files:
        filename = os.path.basename(fasta_file)
        output_tsv = os.path.join(args.out_dir, f"{filename}_blast.tsv")
        tasks.append((fasta_file, args.ref_db, output_tsv, args.blastp_bin, args.threads))

    completed = 0
    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(run_blastp_single, t) for t in tasks]
        for future in as_completed(futures):
            filename, success = future.result()
            completed += 1
            if completed % 20 == 0 or completed == len(fasta_files):
                print(f"Progression : {completed}/{len(fasta_files)} fichiers traités...")

    print(f"\nÉtape 2 terminée : {completed} fichiers TSV générés dans {args.out_dir}")
