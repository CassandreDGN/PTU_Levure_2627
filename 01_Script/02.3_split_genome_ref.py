import os
from pathlib import Path

# Chemins d'entrée et de sortie
input_file = Path("/data/projet1/00_Data/ref_seq_proteique_fasta/orf_trans_all.clean.fasta")
output_dir = Path("/data/projet1/02_Pipeline_Output/2.2_genome_reference")

# Création du dossier de destination
output_dir.mkdir(parents=True, exist_ok=True)

current_file = None
compteur = 0


with open(input_file, 'r') as f_in:
    for line in f_in:
        if line.startswith(">"):
            # Extraction du nom du gène (ex: ">YAL068C description..." -> "YAL068C")
            gene_name = line[1:].strip().split()[0]
            
            # Fermer le fichier précédent s'il y en a un d'ouvert
            if current_file:
                current_file.close()
            
            # Ouvrir un nouveau fichier pour le gène courant
            out_path = output_dir / f"{gene_name}.fasta"
            current_file = open(out_path, 'w')
            current_file.write(line)
            compteur += 1
        else:
            # Écrire la séquence dans le fichier actuellement ouvert
            if current_file:
                current_file.write(line)

# S'assurer de fermer le dernier fichier
if current_file:
    current_file.close()

print(f"Succès ! {compteur} fichiers .fasta ont été générés dans le dossier :")
print(output_dir)