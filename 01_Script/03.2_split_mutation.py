import os
from pathlib import Path

# Définition des chemins d'entrée et de sortie
input_dir = Path("/data/projet1/02_Pipeline_Output/3_Mutations")
output_base_dir = Path("/data/projet1/02_Pipeline_Output/4_Mutations_par_genes")

# Création du dossier principal de sortie s'il n'existe pas
output_base_dir.mkdir(parents=True, exist_ok=True)

# Parcours de tous les fichiers .var dans le dossier d'entrée
for var_file in input_dir.glob("*.var"):
    
    # Extraction du nom de la souche (ex: "AAC.HP1.nuclear_genome" à partir de "AAC.HP1.nuclear_genome.Final...")
    # On coupe la chaîne de caractères au niveau de ".Final"
    strain_name = var_file.name.split(".Final")[0]
    
    # Création du sous-dossier pour cette souche
    strain_out_dir = output_base_dir / strain_name
    strain_out_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Traitement de {var_file.name} -> {strain_name}/")
    
    # Dictionnaire pour stocker les lignes par gène en mémoire (beaucoup plus rapide)
    mutations_by_gene = {}
    
    # Lecture du fichier .var
    with open(var_file, 'r') as f_in:
        for line in f_in:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) >= 2:
                gene_name = parts[0]
                mutation = parts[1] # On isole uniquement la mutation (ex: G432D)
                
                if gene_name not in mutations_by_gene:
                    mutations_by_gene[gene_name] = []
                
                mutations_by_gene[gene_name].append(mutation) # On ajoute que la mutation
    
    # Écriture des fichiers pour chaque gène de cette souche
    for gene, lines in mutations_by_gene.items():
        gene_out_file = strain_out_dir / f"{gene}.var"
        
        # Le mode 'w' écrase le fichier s'il existe déjà, ce qui est idéal si vous relancez le script
        with open(gene_out_file, 'w') as f_out:
            f_out.write('\n'.join(lines) + '\n')

print("\n--- Séparation de tous les fichiers terminée avec succès ! ---")  