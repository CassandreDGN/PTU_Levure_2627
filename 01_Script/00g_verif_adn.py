#!/usr/bin/env python3
"""Que contient-on reellement ? Chercher des donnees ADN (necessaires pour
distinguer faux-sens / synonyme) et les artefacts du scoring."""
import os, glob

print("=== 1. Y a-t-il des donnees ADN / VCF quelque part ? ===")
for motif in ["*.vcf", "*.vcf.gz", "*.gff*", "*.bed", "*cds*.fa", "*.fa", "*.fasta"]:
    r = os.popen(f"find /data/projet1 -name '{motif}' -not -path '*/06_Conda/*' -not -path '*/.git/*' 2>/dev/null | head -6").read().strip()
    if r:
        n = os.popen(f"find /data/projet1 -name '{motif}' -not -path '*/06_Conda/*' -not -path '*/.git/*' 2>/dev/null | wc -l").read().strip()
        print(f"\n[{motif}] -> {n} fichiers")
        print(r[:900])

print("\n=== 2. Contenu d'un fichier 5_Blosum_Scores (artefact du scoring) ===")
f = sorted(glob.glob("/data/projet1/02_Pipeline_Output/5_Blosum_Scores/*"))
if f:
    print("fichier :", os.path.basename(f[0]))
    print(open(f[0]).read()[:900])

print("\n=== 3. 6_Diploid_Genotypes : comment le statut est-il decide ? ===")
f2 = sorted(glob.glob("/data/projet1/02_Pipeline_Output/6_Diploid_Genotypes/*diploid.tsv"))
if f2:
    print("fichier :", os.path.basename(f2[0]))
    print(open(f2[0]).read()[:600])

print("\n=== 4. distribution des phenotypes diploides ===")
r4 = os.popen("cut -f5 /data/projet1/02_Pipeline_Output/6_Diploid_Genotypes/AAC_diploid.tsv 2>/dev/null | sort | uniq -c | sort -rn").read()
print(r4)

print("=== 5. Le script 04.diploid_results.py : sur quoi se base-t-il ? ===")
print(open("/data/projet1/01_Script/04.diploid_results.py").read()[:1800])
