#!/usr/bin/env python3
"""Verification des faits : qui a fait quoi, et quelle est la chronologie reelle."""
import os, glob

print("=== 1. Fichiers dans 04_Resultats_BlastN : proprietaire + date ===")
r = os.popen("ls -la --time-style=long-iso /data/projet1/04_Resultats_BlastN/").read()
print(r)

print("=== 2. Dates cles du projet (auteurs differencies) ===")
r2 = os.popen(
    "for f in /data/projet1/01_Script/*.py /data/projet1/01_Script/*.sh; do "
    "stat -c '%y  %U  %n' $f; done | sort").read()
print(r2)

print("=== 3. Les resultats de l'etape 3 : qui les a produits ? ===")
for d in ["3_Mutations", "5_Blosum_Scores", "6_Diploid_Genotypes", "7_Synthesis_and_Plots"]:
    p = f"/data/projet1/02_Pipeline_Output/{d}"
    if os.path.isdir(p):
        fs = sorted(glob.glob(p + "/*"))
        if fs:
            print(f"{d}: {len(fs)} fichiers, premier = {os.path.basename(fs[0])}, proprietaire = {os.popen(f'stat -c %U {fs[0]}').read().strip()}, date = {os.popen(f'stat -c %y {fs[0]}').read().strip()[:16]}")

print("\n=== 4. Le script 03.1 (extraction alignee) : auteur ===")
print(os.popen("stat -c '%y  %U  %n' /data/projet1/01_Script/03.1_extract_mutations.py").read())

print("=== 5. Les .var de 3_Mutations : proprietaire + date ===")
print(os.popen("ls -la --time-style=long-iso /data/projet1/02_Pipeline_Output/3_Mutations/ | head -5").read())
print(os.popen("ls -la --time-style=long-iso /data/projet1/02_Pipeline_Output/3_Mutations/ | tail -3").read())
