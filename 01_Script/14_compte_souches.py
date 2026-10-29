#!/usr/bin/env python3
"""Compte les souches de 3 lettres vs les cas particuliers (SACE_*, XTRA_*, noms sans HP)."""
import csv
from collections import defaultdict

SRC = "/data/projet1/02_Pipeline_Output/8_Stats_Mutations/stats_mutations_par_souche.tsv"

donnees = defaultdict(dict)
with open(SRC) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        donnees[row["souche"]][row["haplotype"]] = row["total_variants"]

trois = [s for s in donnees if len(s) == 3 and not s.startswith(("SACE_", "XTRA_"))]
sace = [s for s in donnees if s.startswith("SACE_")]
xtra = [s for s in donnees if s.startswith("XTRA_")]
sans_hp = [s for s in donnees if "NA" in donnees[s] or len(donnees[s]) == 0]

print(f"Souches 3 lettres (cas general) : {len(trois)}")
print(f"Souches SACE_* : {len(sace)}")
print(f"Souches XTRA_* : {len(xtra)}")
print(f"Souches sans haplotype distinct : {len(sans_hp)} -> {sorted(sans_hp)[:10]}")
print(f"TOTAL : {len(donnees)}")

# souches avec les 2 haplotypes
deux = [s for s in donnees if len(donnees[s]) == 2]
print(f"\nSouches avec HP1 + HP2 : {len(deux)}")

# et les fichiers sans .HP1/.HP2 dans le nom ?
import glob, os
fichiers = sorted(os.path.basename(f) for f in glob.glob("/data/projet1/02_Pipeline_Output/3_Mutations/*.var"))
anomalie = [f for f in fichiers if ".HP1." not in f and ".HP2." not in f]
print(f"\nFichiers .var sans .HP1/.HP2 dans le nom : {len(anomalie)}")
for f in anomalie:
    print("   ", f)
