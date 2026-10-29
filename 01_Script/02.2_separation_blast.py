import os
from pathlib import Path

# Définition des chemins d'entrée et de sortie
input_dir = Path("/data/projet1/02_Pipeline_Output/2_Blast")
output_base_dir = Path("/data/projet1/02_Pipeline_Output/2.1_blast_par_genes")

# Création du dossier principal de sortie s'il n'existe pas
output_base_dir.mkdir(parents=True, exist_ok=True)

# Parcours de tous les fichiers TSV dans le dossier d'entrée
for blast_file in input_dir.glob("*_blast.tsv"):
    
    # Extraction du nom de la souche (ex: "AAC.HP1.nuclear_genome")
    strain_name = blast_file.name.split(".Final")[0]
    
    # Création du sous-dossier pour cette souche
    strain_out_dir = output_base_dir / strain_name
    strain_out_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Traitement de {blast_file.name} -> {strain_name}/")
    
    # Dictionnaire pour stocker les lignes par gène
    blast_by_gene = {}
    
    # Lecture du fichier BLAST
    with open(blast_file, 'r') as f_in:
        for line in f_in:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) >= 2:
                # Dans un fichier BLAST au format tabulaire, le gène (sujet) est la 2ème colonne (index 1)
                gene_name = parts[1]
                
                if gene_name not in blast_by_gene:
                    blast_by_gene[gene_name] = []
                
                blast_by_gene[gene_name].append(line)
    
    # Écriture des fichiers pour chaque gène
    for gene, lines in blast_by_gene.items():
        # J'ai nommé les fichiers avec l'extension _blast.tsv pour garder leur nature explicite.
        # Si vous préférez strictement ".fasta" ou ".tsv", modifiez simplement la ligne ci-dessous.
        gene_out_file = strain_out_dir / f"{gene}_blast.tsv"
        
        with open(gene_out_file, 'w') as f_out:
            f_out.write('\n'.join(lines) + '\n')

print("\n--- Séparation de tous les fichiers BLAST terminée avec succès ! ---")