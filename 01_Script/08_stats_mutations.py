#!/usr/bin/env python3
"""Analyse statistique de l'etape 3 (mutations) sur les fichiers .var de toutes les souches.
Produit un tableau par souche + un resume global. Ne depend que de la stdlib."""
import os, re, glob, csv, json
from collections import Counter, defaultdict

VAR_DIR = "/data/projet1/02_Pipeline_Output/3_Mutations"
OUT_DIR = "/data/projet1/02_Pipeline_Output/8_Stats_Mutations"
os.makedirs(OUT_DIR, exist_ok=True)

# formes attendues : "<AA><pos><AA>" substitution, "<AA><pos>del", "<p1>_<p2>ins<AA>"
RE_SUB = re.compile(r"^([A-Z*])(\d+)([A-Z*])$")
RE_DEL = re.compile(r"^([A-Z*])(\d+)del$")
RE_INS = re.compile(r"^(\d+)_(\d+)ins([A-Z*]+)$")

var_files = sorted(glob.glob(os.path.join(VAR_DIR, "*.var")))
print(f"Fichiers .var trouves : {len(var_files)}")

per_strain = []          # lignes du tableau resume
tot_type = Counter()     # types de variants (global)
tot_aa = Counter()       # substitutions AA->AA (global)
gene_hits = Counter()    # gènes les plus mutés
tot_lines = 0
genes_total = 0

for vf in var_files:
    base = os.path.basename(vf).replace(".nuclear_genome.Final.cds.prot.fasta.var", "")
    strain = base.split(".")[0] if "." in base else base
    hap = base.split(".")[1] if "." in base else "NA"

    n_sub = n_del = n_ins = 0
    genes = set()
    aa = Counter()
    with open(vf) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            gene, var = parts[0], parts[1]
            tot_lines += 1
            genes.add(gene)
            m = RE_SUB.match(var)
            if m:
                n_sub += 1
                tot_type["substitution"] += 1
                aa[(m.group(1), m.group(3))] += 1
                tot_aa[(m.group(1), m.group(3))] += 1
                continue
            if RE_DEL.match(var):
                n_del += 1
                tot_type["deletion"] += 1
                continue
            if RE_INS.match(var):
                n_ins += 1
                tot_type["insertion"] += 1
                continue
            tot_type["autre"] += 1

    n_tot = n_sub + n_del + n_ins
    per_strain.append({
        "souche": strain, "haplotype": hap, "total_variants": n_tot,
        "genes_mutés": len(genes), "substitutions": n_sub,
        "deletions": n_del, "insertions": n_ins,
        "variants_par_gene": round(n_tot / len(genes), 2) if genes else 0,
    })
    genes_total += len(genes)
    print(f"  {base:22s} total={n_tot:7d} sub={n_sub:7d} del={n_del:5d} ins={n_ins:5d} genes={len(genes)}", flush=True)

# --- tableau par souche ---
per_strain.sort(key=lambda d: -d["total_variants"])
with open(os.path.join(OUT_DIR, "stats_mutations_par_souche.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(per_strain[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(per_strain)

# --- resume global ---
n_files = len(per_strain)
totals = [d["total_variants"] for d in per_strain]
summary = {
    "n_fichiers_var": n_files,
    "n_souches": len({d["souche"] for d in per_strain}),
    "total_variants": sum(totals),
    "moyenne_variants_par_fichier": round(sum(totals) / n_files, 1),
    "mediane_variants_par_fichier": sorted(totals)[n_files // 2],
    "min": min(totals), "max": max(totals),
    "types": dict(tot_type),
    "pct_substitution": round(100 * tot_type["substitution"] / sum(tot_type.values()), 1),
    "pct_deletion": round(100 * tot_type["deletion"] / sum(tot_type.values()), 1),
    "pct_insertion": round(100 * tot_type["insertion"] / sum(tot_type.values()), 1),
}
with open(os.path.join(OUT_DIR, "resume_mutations_global.json"), "w") as fh:
    json.dump(summary, fh, indent=2, ensure_ascii=False)

print("\n=== RESUME GLOBAL ===")
print(json.dumps(summary, indent=2, ensure_ascii=False))

print("\n=== TOP 20 SUBSTITUTIONS AA ===")
for (a, b), c in tot_aa.most_common(20):
    print(f"  {a} -> {b} : {c}")

print("\n=== TOP 15 SOUCHES (nb variants) ===")
for d in per_strain[:15]:
    print(f"  {d['souche']:10s} {d['haplotype']:4s} total={d['total_variants']:7d} sub={d['substitutions']:7d} "
          f"del={d['deletions']:5d} ins={d['insertions']:5d} genes={d['genes_mutés']}")

print("\n=== TOP 15 SOUCHES (nb genes mutés) ===")
for d in sorted(per_strain, key=lambda x: -x["genes_mutés"])[:15]:
    print(f"  {d['souche']:10s} genes_mutés={d['genes_mutés']:5d} total_variants={d['total_variants']:7d}")

print("\nSorties dans", OUT_DIR)
