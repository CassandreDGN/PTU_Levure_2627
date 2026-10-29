#!/usr/bin/env python3
"""QC global de l'etape 1 (traduction) : lit tous les qc_summary.tsv et compte
les CDS rejetes par motif (not_multiple_3, no_atg, invalid_codon, premature_stop,
too_short, empty). Produit la table qui manque au cahier de labo."""
import os, glob, csv
from collections import defaultdict

ROOT = "/data/projet1/02_Pipeline_Output"
OUT = "/data/projet1/02_Pipeline_Output/8_Stats_Mutations/qc_traduction_global.tsv"

MOTIFS = ["total_cds", "valid_prot", "not_multiple_3", "no_atg",
          "invalid_codon", "premature_stop", "too_short", "empty"]

tot = defaultdict(int)
rows = []
files = sorted(glob.glob(os.path.join(ROOT, "**", "qc_summary.tsv"), recursive=True))
for f in files:
    with open(f) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            rows.append((os.path.basename(os.path.dirname(f)), row["souche"], row))
            for m in MOTIFS:
                tot[m] += int(row.get(m) or 0)

with open(OUT, "w", newline="") as fh:
    fh.write("source\tsouche\t" + "\t".join(MOTIFS) + "\n")
    for src, s, row in rows:
        fh.write(f"{src}\t{s}\t" + "\t".join(str(row.get(m) or 0) for m in MOTIFS) + "\n")

print(f"Fichiers qc_summary.tsv trouves : {len(files)}")
print(f"Lignes (souches) : {len(rows)}\n")
print("=== TOTAUX QC TRADUCTION ===")
for m in MOTIFS:
    print(f"  {m:18s} {tot[m]}")
print(f"\nTaux de rejet total : {100*(tot['total_cds']-tot['valid_prot'])/max(tot['total_cds'],1):.2f} %")
print("Sortie :", OUT)
