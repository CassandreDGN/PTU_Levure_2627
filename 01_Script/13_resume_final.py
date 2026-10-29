#!/usr/bin/env python3
"""Chiffres manquants pour le resume : moyenne par SOUCHE diploide (2 haplotypes
confondus) et detail par haplotype. Ecrit dans 8_Stats_Mutations."""
import csv, json
from collections import defaultdict

SRC = "/data/projet1/02_Pipeline_Output/8_Stats_Mutations/stats_mutations_par_souche.tsv"
OUT = "/data/projet1/02_Pipeline_Output/8_Stats_Mutations/resume_final_mutations.md"

# souche -> {hp: ligne}
donnees = defaultdict(dict)
with open(SRC) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        s = row["souche"]
        if s.startswith("SACE_") or s.startswith("XTRA_"):
            continue
        donnees[s][row["haplotype"]] = {
            "total": int(row["total_variants"]),
            "genes": int(row["genes_mutés"]),
            "sub": int(row["substitutions"]),
            "del": int(row["deletions"]),
            "ins": int(row["insertions"]),
        }

par_souche = []
for s, hp in donnees.items():
    if len(hp) != 2:
        continue
    tot = sum(v["total"] for v in hp.values())
    g = max(v["genes"] for v in hp.values())
    par_souche.append((s, tot, g, hp["HP1"]["total"], hp["HP2"]["total"]))

totaux = [t for _, t, _, _, _ in par_souche]
n = len(par_souche)
moy = sum(totaux) / n
tries = sorted(totaux)

lignes = []
lignes.append("# Etape 3 — Extraction et caractérisation des mutations\n")
lignes.append(f"Periode : 19-21 octobre 2026. Auteur : Beray Sonmez.\n")
lignes.append("## 1. Ce qui a ete fait\n")
lignes.append(
    "Apres l'alignement de chaque CDS de souche sur son gene de reference S288C "
    "(etape 2, DIAMOND/BLASTp), les differences ont ete extraites acide amine par "
    "acide amine et ecrites au format de variations. Trois classes de variants sont "
    "produites : substitutions, deletions et insertions.\n")
lignes.append("## 2. Resultats\n")
lignes.append(f"- Fichiers de variations produits : **797** (397 souches x 2 haplotypes)")
lignes.append(f"- Variants totaux : **24 809 068**")
lignes.append(f"- Substitutions : **16 242 651** (65,5 %)")
lignes.append(f"- Deletions : **7 063 762** (28,5 %)")
lignes.append(f"- Insertions : **1 502 655** (6,1 %)\n")
lignes.append(f"### Moyenne par souche diploide\n")
lignes.append(f"- **{moy:.0f} variants par souche** (les 2 haplotypes reunis)")
lignes.append(f"- Mediane : {tries[n//2]} variants")
lignes.append(f"- Etendue : {tries[0]} a {tries[-1]} variants")
lignes.append(f"- Nombre de souches avec 2 haplotypes : {n}\n")
lignes.append("### Variants par type, en moyenne par souche diploide\n")
lignes.append(f"- {16_242_651/n:.0f} substitutions")
lignes.append(f"- {7_063_762/n:.0f} deletions")
lignes.append(f"- {1_502_655/n:.0f} insertions\n")

with open(OUT, "w") as fh:
    fh.write("\n".join(lignes) + "\n")

print(f"n souches diploides = {n}")
print(f"moyenne par souche diploide = {moy:.0f}")
print(f"mediane = {tries[n//2]}")
print(f"min = {tries[0]}  max = {tries[-1]}")
print("\n--- top 10 souches ---")
for s, t, g, h1, h2 in sorted(par_souche, key=lambda x: -x[1])[:10]:
    print(f"  {s:6s} total={t:7d} (HP1={h1:6d} HP2={h2:6d}) genes={g}")
print("\nSortie :", OUT)
