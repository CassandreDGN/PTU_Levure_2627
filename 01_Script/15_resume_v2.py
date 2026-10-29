#!/usr/bin/env python3
"""Resume final corrige : souches 3 lettres avec 2 haplotypes uniquement (358),
et mention des cas particuliers. Ecrase resume_final_mutations.md."""
import csv
from collections import defaultdict

SRC = "/data/projet1/02_Pipeline_Output/8_Stats_Mutations/stats_mutations_par_souche.tsv"
OUT = "/data/projet1/02_Pipeline_Output/8_Stats_Mutations/resume_final_mutations.md"

donnees = defaultdict(dict)
with open(SRC) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        s = row["souche"]
        if len(s) == 3 and not s.startswith(("SACE_", "XTRA_")):
            donnees[s][row["haplotype"]] = int(row["total_variants"])

complet = {s: v for s, v in donnees.items() if len(v) == 2}
totaux = [sum(v.values()) for v in complet.values()]
n = len(complet)
tries = sorted(totaux)

L = []
L.append("# Etape 3 — Extraction et caracterisation des mutations")
L.append("")
L.append("Periode : 19-21 octobre 2026. Auteur : Beray Sonmez.")
L.append("")
L.append("## 1. Ce qui a ete fait")
L.append("")
L.append("Apres l'alignement de chaque CDS de souche sur son gene de reference S288C")
L.append("(etape 2, DIAMOND/BLASTp), les differences ont ete extraites acide amine par")
L.append("acide amine et ecrites au format de variations. Trois classes de variants sont")
L.append("produites : substitutions, deletions, insertions.")
L.append("")
L.append("## 2. Resultats globaux")
L.append("")
L.append("- Fichiers de variations produits : **797** (397 souches x 2 haplotypes)")
L.append("- Variants totaux : **24 809 068**")
L.append("- Substitutions : **16 242 651** (65,5 %)")
L.append("- Deletions : **7 063 762** (28,5 %)")
L.append("- Insertions : **1 502 655** (6,1 %)")
L.append("")
L.append("## 3. Moyenne par souche diploide")
L.append("")
L.append(f"Base : {n} souches de 3 lettres possedant les 2 haplotypes (HP1 + HP2).")
L.append(f"Les 17 souches SACE_* et 19 souches XTRA_* sont exclues de ce calcul, ainsi que")
L.append("les 3 souches sans fichier haplotype distinct (AHP, BHP, CHP).")
L.append("")
L.append(f"- **{sum(totaux)/n:.0f} variants par souche diploide** (2 haplotypes reunis)")
L.append(f"- Mediane : {tries[n//2]} variants")
L.append(f"- Etendue : {tries[0]} a {tries[-1]} variants")
L.append("")
L.append("En moyenne par souche diploide et par type :")
L.append(f"- {16_242_651/n:.0f} substitutions")
L.append(f"- {7_063_762/n:.0f} deletions")
L.append(f"- {1_502_655/n:.0f} insertions")
L.append("")
L.append("## 4. Souches les plus variables")
L.append("")
L.append("| Souche | Total | HP1 | HP2 |")
L.append("|---|---|---|---|")
for s, v in sorted(complet.items(), key=lambda x: -sum(x[1].values()))[:8]:
    L.append(f"| {s} | {sum(v.values())} | {v['HP1']} | {v['HP2']} |")
L.append("")
L.append("## 5. Controle qualite")
L.append("")
L.append("- Variants vers un codon STOP : **0**")
L.append("- Variants a position invalide (au-dela de la proteine de reference) : **0**")
L.append("- Insertions de plus d'un acide amine : **0**")
L.append("")
L.append("Ces trois compteurs a zero montrent que les variants extraits sont reels et")
L.append("que la traduction en amont a bien ecarte les sequences a codon stop premature.")
L.append("")
L.append("## 6. Bilan de l'annotation (rattachement CDS -> gene de reference)")
L.append("")
L.append("- CDS traduits : 4 719 625")
L.append("- CDS rattaches a un gene S288C : 4 652 953 (**99,3 %**)")
L.append("- Moyenne par souche x haplotype : 5 838 genes rattaches")
L.append("- Genes de reference dans lesquels des variants ont ete trouves : 5 898")
L.append("")
L.append("## 7. Fichiers produits")
L.append("")
L.append("- `stats_mutations_par_souche.tsv` — 797 lignes (souche, haplotype, totaux)")
L.append("- `stats_mutations_par_gene.tsv` — 5 898 genes classes par nombre de variants")
L.append("- `annotation_par_souche.tsv` — genes rattaches par souche")
L.append("- `resume_mutations_global.json` — totaux et pourcentages")
L.append("- `sanity_variants.tsv` — compteurs de controle qualite")
L.append("- `qc_traduction_global.tsv` — QC de l'etape 1")

with open(OUT, "w") as fh:
    fh.write("\n".join(L) + "\n")

print(f"n = {n} souches 3 lettres avec 2 haplotypes")
print(f"moyenne = {sum(totaux)/n:.0f}")
print("->", OUT)
