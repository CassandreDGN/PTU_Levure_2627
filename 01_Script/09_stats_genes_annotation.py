#!/usr/bin/env python3
"""Etape 3 bis : analyse par gène (lesquels sont les plus touches) +
distribution des tailles + bilan d'annotation (genes de souche rattaches a la reference).
Ecrit dans /data/projet1/02_Pipeline_Output/8_Stats_Mutations."""
import os, re, glob, csv, json
from collections import Counter, defaultdict

VAR_DIR = "/data/projet1/02_Pipeline_Output/3_Mutations"
BLAST_DIR = "/data/projet1/02_Pipeline_Output/2_Blast"
OUT_DIR = "/data/projet1/02_Pipeline_Output/8_Stats_Mutations"
os.makedirs(OUT_DIR, exist_ok=True)

RE_SUB = re.compile(r"^([A-Z*])(\d+)([A-Z*])$")
RE_DEL = re.compile(r"^([A-Z*])(\d+)del$")
RE_INS = re.compile(r"^(\d+)_(\d+)ins([A-Z*]+)$")

# ---------- 1. mutations par gène (toutes souches confondues) ----------
gene_var = Counter()
gene_strains = defaultdict(set)
for vf in sorted(glob.glob(os.path.join(VAR_DIR, "*.var"))):
    base = os.path.basename(vf).replace(".nuclear_genome.Final.cds.prot.fasta.var", "")
    strain = base.split(".")[0] if "." in base else base
    with open(vf) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            gene_var[p[0]] += 1
            gene_strains[p[0]].add(strain)

with open(os.path.join(OUT_DIR, "stats_mutations_par_gene.tsv"), "w", newline="") as fh:
    fh.write("gene_reference\tn_variants\tn_souches_muté\n")
    for g, c in gene_var.most_common():
        fh.write(f"{g}\t{c}\t{len(gene_strains[g])}\n")

print(f"Genes differents porteurs de variants : {len(gene_var)}")
print("\n=== TOP 20 GENES (nb variants) ===")
for g, c in gene_var.most_common(20):
    print(f"  {g:16s} {c:8d} variants, {len(gene_strains[g]):3d} souches")

# ---------- 2. bilan d'annotation : CDS de souche -> gene de reference ----------
tsv = sorted(glob.glob(os.path.join(BLAST_DIR, "*_blast.tsv")))
n_cds_total = 0
n_cds_aligned = 0
per_strain_cds = []
for f in tsv:
    base = os.path.basename(f).replace("_blast.tsv", "")
    ids = set()
    with open(f) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            ids.add(p[0])
    per_strain_cds.append((base, len(ids)))
    n_cds_aligned += len(ids)

# nombre de CDS dans les FASTA de souche traduits
prot_files = sorted(glob.glob("/data/projet1/02_Pipeline_Output/1_Translated_Proteins/*.prot.fasta"))
for pf in prot_files:
    with open(pf) as fh:
        n_cds_total += sum(1 for l in fh if l.startswith(">"))

with open(os.path.join(OUT_DIR, "annotation_par_souche.tsv"), "w", newline="") as fh:
    fh.write("souche_haplotype\tn_genes_alignes\n")
    for b, c in sorted(per_strain_cds):
        fh.write(f"{b}\t{c}\n")

print(f"\n=== BILAN ANNOTATION ===")
print(f"Fichiers tsv (souches x haplotypes) : {len(tsv)}")
print(f"CDS traduits (total, tous fichiers) : {n_cds_total}")
print(f"Genes de souche alignes sur la reference : {n_cds_aligned}")
print(f"Genes differents dans la reference : 6722")
print(f"Taux de rattachement : {100*n_cds_aligned/n_cds_total:.1f}%")

moy = n_cds_aligned / max(len(tsv), 1)
print(f"Moyenne par souche x haplotype : {moy:.0f} genes rattaches")
print(f"\nSorties dans {OUT_DIR}")
