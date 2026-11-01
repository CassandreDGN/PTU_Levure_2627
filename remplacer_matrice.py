#!/usr/bin/env python3
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

# Chemin vers tes résultats Swiss-Prot
TARGET_DIR = Path("/data/projet1/02_Pipeline_Output/4_SIFT_SWISSPROT_Results")

def clean_sift_file(file_path):
    try:
        # 1. On lit seulement le début du fichier (ultra-rapide)
        with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
            header = f.read(150)
        
        # 2. Si c'est une matrice d'erreur, on écrase par "matrice"
        if "UNK_ID; MATRIX" in header or "MATRIX" in header and "UNK_ID" in header:
            with open(file_path, "w", encoding="utf-8") as f:
                f.write("matrice\n")
            return 1
    except Exception:
        pass
    return 0

if __name__ == "__main__":
    print("=== Nettoyage des matrices SIFT (Swiss-Prot) ===")
    
    # Récupération de tous les fichiers .SIFTprediction
    sift_files = list(TARGET_DIR.rglob("*.SIFTprediction"))
    print(f"Nombre total de fichiers SIFT à vérifier : {len(sift_files)}")

    # Exécution sur x cœurs en parallèle
    with ProcessPoolExecutor(max_workers=5) as executor:
        results = executor.map(clean_sift_file, sift_files, chunksize=1000)
        modified = sum(results)

    print(f"Terminé ! {modified} fichiers ont été remplacés par 'matrice'.")