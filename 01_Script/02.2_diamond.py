#!/usr/bin/env python3
import os
import glob
import argparse
import subprocess
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed

def run_diamond_single(args_tuple):
    strain_fasta, ref_db, output_tsv, diamond_bin, threads = args_tuple
    filename = os.path.basename(strain_fasta)
    
    cmd = [
    diamond_bin, "blastp",
    "--query", strain_fasta,
    "--db", ref_db,
    "--out", output_tsv,
    "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "qseq", "sseq", "sstart",
    "--max-target-seqs", "1",
    "--evalue", "1e-5",
    "--masking", "0",  # <--- Désactive le masquage par X
    "--threads", str(threads)
]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        return filename, True, None
    except subprocess.CalledProcessError as e:
        return filename, False, e.stderr
    except Exception as e:
        return filename, False, str(e)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Exécute DIAMOND blastp en parallèle.")
    
    parser.add_argument("-i", "--data_dir", default="/data/projet1/02_Pipeline_Output/1_Translated_Proteins")
    parser.add_argument("-db", "--ref_db", default="/data/projet1/00_Data/db_diamond/Ref_S288C_Prot.dmnd")
    parser.add_argument("-o", "--out_dir", default="/data/projet1/02_Pipeline_Output/2_Blast")
    parser.add_argument("-b", "--diamond_bin", default=shutil.which("diamond") or "/home/degen/.conda/envs/blastp_env/bin/diamond")
    parser.add_argument("-j", "--jobs", type=int, default=4)
    parser.add_argument("-t", "--threads", type=int, default=1)

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # Vérification de la base DIAMOND
    if not os.path.exists(args.ref_db):
        print(f"Erreur : La base DIAMOND est introuvable à l'adresse {args.ref_db}", flush=True)
        print("Générez-la d'abord avec : diamond makedb --in <ref.fasta> -d /data/projet1/00_Data/db_diamond/Ref_S288C_Prot", flush=True)
        exit(1)

    fasta_files = sorted(glob.glob(os.path.join(args.data_dir, "*.prot.fasta")))

    if not fasta_files:
        print(f"Erreur : Aucun fichier .prot.fasta trouvé dans {args.data_dir}", flush=True)
        exit(1)

    print(f"Lancement de DIAMOND blastp sur {len(fasta_files)} fichiers .prot.fasta...", flush=True)
    print(f"Base de données : {args.ref_db}", flush=True)
    print(f"Exécutable : {args.diamond_bin}\n", flush=True)

    tasks = []
    for fasta_file in fasta_files:
        filename = os.path.basename(fasta_file)
        output_tsv = os.path.join(args.out_dir, f"{filename}_blast.tsv")
        tasks.append((fasta_file, args.ref_db, output_tsv, args.diamond_bin, args.threads))

    completed = 0
    errors = 0
    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(run_diamond_single, t) for t in tasks]
        for future in as_completed(futures):
            filename, success, err_msg = future.result()
            if success:
                completed += 1
            else:
                errors += 1
                print(f"ÉCHEC pour {filename} : {err_msg}", flush=True)
            
            if (completed + errors) % 20 == 0 or (completed + errors) == len(fasta_files):
                print(f"Progression : {completed + errors}/{len(fasta_files)} traités ({completed} succès, {errors} échecs)...", flush=True)

    print(f"\nÉtape 2 terminée : {completed} fichiers TSV générés avec succès dans {args.out_dir}", flush=True)