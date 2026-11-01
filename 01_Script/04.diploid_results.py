import re
from pathlib import Path
from collections import defaultdict

# CORRECTION : Pointer vers le dossier où 03.5 a créé les fichiers summary
IN_DIR = Path("/data/projet1/02_Pipeline_Output/3_SIFT_Mutations")
OUT_DIR = Path("/data/projet1/02_Pipeline_Output/6_Diploid_Genotypes")
OUT_DIR.mkdir(parents=True, exist_ok=True)

print("=== ÉTAPE 4 : ÉVALUATION DES GÉNOTYPES DIPLOÏDES ===")

strain_copies = defaultdict(dict)

# Diagnostic des dossiers
all_folders = [f for f in IN_DIR.iterdir() if f.is_dir()]
print(f"Dossiers détectés dans {IN_DIR.name} : {len(all_folders)}")

for folder in all_folders:
    match = re.search(r'^(.*?)\.HP([12])(?:\..*)?$', folder.name)
    if match:
        strain_name = match.group(1)
        copy_num = match.group(2)
        summary_file = folder / "gene_summary.tsv"
        
        if summary_file.exists():
            strain_copies[strain_name][copy_num] = summary_file
        # On ne met pas de warning ici pour les dossiers en attente (sinon tu auras 796 avertissements)

# Nettoyage des souches incomplètes
complete_strains = {name: copies for name, copies in strain_copies.items() if '1' in copies and '2' in copies}
print(f"-> {len(complete_strains)} souches diploïdes complètes identifiées.")

global_matrix = defaultdict(dict)
all_genes = set()

for strain_name, copies in sorted(complete_strains.items()):
    status_copy1, status_copy2 = {}, {}

    # Lecture HP1
    with open(copies['1']) as f:
        f.readline()  # Saut de l'en-tête
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 7:
                # parts[-1] extrait 'is_functional' (YES/NO/UNKNOWN)
                status_copy1[parts[0]] = 1 if parts[-1] == 'YES' else 0

    # Lecture HP2
    with open(copies['2']) as f:
        f.readline()  # Saut de l'en-tête
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 7:
                status_copy2[parts[0]] = 1 if parts[-1] == 'YES' else 0

    strain_genes = set(status_copy1.keys()).union(set(status_copy2.keys()))
    all_genes.update(strain_genes)

    out_strain_file = OUT_DIR / f"{strain_name}_diploid.tsv"
    
    functional_homo = 0
    functional_hetero = 0
    total_loss = 0

    with open(out_strain_file, "w") as fout:
        fout.write("protein\thp1_status\thp2_status\tdiploid_score\tphenotype\n")
        
        for prot_id in sorted(strain_genes):
            c1 = status_copy1.get(prot_id, 0) # 0 par défaut si manquant
            c2 = status_copy2.get(prot_id, 0)
            
            diploid_score = c1 + c2
            
            if diploid_score == 2:
                pheno = "FUNCTIONAL"
                functional_homo += 1   # Tu peux renommer la variable functional_homo en count_functional si tu veux
            elif diploid_score == 1:
                pheno = "MONOALLELIC_LOF"
                functional_hetero += 1 # Tu peux renommer en count_monoallelic
            elif diploid_score == 0:
                pheno = "BIALLELIC_LOF"
                total_loss += 1

            fout.write(f"{prot_id}\t{c1}\t{c2}\t{diploid_score}\t{pheno}\n")
            global_matrix[prot_id][strain_name] = diploid_score

    print(f"  [OK] Souche {strain_name} -> Fonctionnel (2): {functional_homo} | Perte Monoallélique (1): {functional_hetero} | LoF Biallélique (0): {total_loss}")

# Export de la matrice globale gènes x souches
matrix_file = OUT_DIR / "diploid_matrix_all_strains.tsv"
if all_genes:
    sorted_strains = sorted(next(iter(global_matrix.values())).keys())
    with open(matrix_file, "w") as fmat:
        fmat.write("protein\t" + "\t".join(sorted_strains) + "\n")
        for prot_id in sorted(all_genes):
            # Le score par défaut dans la matrice devient 0 (si la souche n'a aucune donnée pour ce gène)
            row_scores = [str(global_matrix[prot_id].get(s, 0)) for s in sorted_strains]
            fmat.write(f"{prot_id}\t" + "\t".join(row_scores) + "\n")

    print(f"\n Matrice diploïde globale générée : {matrix_file.name}")
else:
    print("\n Aucune donnée à exporter. Vérifiez les fichiers.")