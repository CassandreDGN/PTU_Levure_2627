    #!/usr/bin/env python3
    import os
    import glob
    import argparse
    import subprocess

    def run_blastp_souche(strain_fasta, ref_db, output_tsv):
        cmd = [
            "blastp",
            "-query", strain_fasta,
            "-db", ref_db,
            "-out", output_tsv,
            "-outfmt", "6 qseqid sseqid pident length evalue bitscore qseq sseq",
            "-max_target_seqs", "1",
            "-evalue", "1e-5"
        ]
        # qseqid   = id de la prot de la souche
        # sseqid   = id du gene de la ref
        # pident   = pourcentage d'identité
        # length   = longueur de l'alignement
        # evalue   = e-value (plus proche de 0, plus l'alignement est significatif)
        # bitscore = score de l'alignement
        # qseq     = sequence d'aa de la souche
        # sseq     = sequence d'aa de la ref
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=None)
            return True
        except subprocess.CalledProcessError as e:
            print(f"Erreur lors du BLAST pour {os.path.basename(strain_fasta)}: {e}")
            return False

    if __name__ == "__main__":
        # Configuration des arguments pour le style
        parser = argparse.ArgumentParser(
            description="Exécute le BLASTp pour chaque souche contre la référence S288C."
        )
        
        parser.add_argument(
            "-i", "--data_dir",
            default="/data/projet1/00_Data/cds_fasta",
            help="Dossier contenant les fichiers FASTA des souches (défaut: %(default)s)"
        )
        parser.add_argument(
            "-db", "--ref_db",
            default="/data/projet1/00_Data/reference_S288C/orf_trans_all",
            help="Chemin vers la base de données BLAST de référence (défaut: %(default)s)"
        )
        parser.add_argument(
            "-o", "--out_dir",
            default="/data/projet1/02_Pipeline_Output/02_Blast",
            help="Dossier de sortie pour les fichiers TSV (défaut: %(default)s)"
        )

        args = parser.parse_args()

        data_dir = args.data_dir
        ref_db = args.ref_db
        out_dir = args.out_dir
        
        os.makedirs(out_dir, exist_ok=True)

        fasta_files = glob.glob(os.path.join(data_dir, "*.fasta"))

        for fasta_file in fasta_files:
            filename = os.path.basename(fasta_file)
            strain_name = os.path.splitext(filename)[0]
            
            output_tsv = os.path.join(out_dir, f"{strain_name}_blast.tsv")
            
            print(f"BLASTp en cours pour la souche : {strain_name}...")
            run_blastp_souche(fasta_file, ref_db, output_tsv)
            
        print(f"\nÉtape 2 terminée : Fichiers TSV générés dans {out_dir}")