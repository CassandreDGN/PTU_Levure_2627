#!/usr/bin/env python3
"""Ce qui decide du statut dans 04.diploid_results.py + contenu des dossiers de scoring."""
import os, glob

print("=== 1. 04.diploid_results.py : logique de decision ===")
print(open("/data/projet1/01_Script/04.diploid_results.py").read())

print("=== 2. Contenu d'un dossier 5_Blosum_Scores ===")
d = "/data/projet1/02_Pipeline_Output/5_Blosum_Scores"
sous = sorted(glob.glob(d + "/*"))
print("nombre d'entrees :", len(sous))
print("premiere entree :", sous[0])
if os.path.isdir(sous[0]):
    for f in sorted(glob.glob(sous[0] + "/*"))[:6]:
        print("   ", os.path.basename(f), os.path.getsize(f), "octets")
    # lire le premier fichier de resultat
    fs = [f for f in sorted(glob.glob(sous[0] + "/*")) if os.path.isfile(f)]
    if fs:
        print("--- contenu de", os.path.basename(fs[0]), "---")
        print(open(fs[0]).read()[:1200])

print("\n=== 3. phenotypes dans 6_Diploid_Genotypes ===")
r = os.popen("cut -f5 /data/projet1/02_Pipeline_Output/6_Diploid_Genotypes/AAC_diploid.tsv 2>/dev/null | sort | uniq -c | sort -rn").read()
print(r)
print("--- extrait ---")
print(os.popen("head -6 /data/projet1/02_Pipeline_Output/6_Diploid_Genotypes/AAC_diploid.tsv").read())

print("=== 4. Y a-t-il un VCF ou un fichier de variants ADN ? ===")
print(os.popen("find /data/projet1 -iname '*vcf*' -o -iname '*snp*' -o -iname '*variant*' 2>/dev/null | grep -v 06_Conda | grep -v .git | head -10").read())
print("(vide = aucun)")

print("=== 5. Les CDS de reference S288C : ADN disponible ? ===")
print(os.popen("head -2 /data/projet1/00_Data/reference_S288C/orf_coding_all.fasta").read()[:300])
print(os.popen("grep -c '^>' /data/projet1/00_Data/reference_S288C/orf_coding_all.fasta").read())
