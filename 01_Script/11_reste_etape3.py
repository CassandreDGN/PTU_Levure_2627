#!/usr/bin/env python3
"""Ce qui reste de l'etape 3 : les variants auxquels on n'a PAS encore applique de
prediction de fonctionnalite (truncations, codons stop premature) + la liste des
protéines de souche non encore traitees. Ecrit dans 02_Pipeline_Output/8_Stats_Mutations.
"""
import os, re, glob, csv, json

REF = "/data/projet1/00_Data/ref_seq_proteique_fasta/orf_trans_all.fasta"
PROT_DIR = "/data/projet1/02_Pipeline_Output/1_Translated_Proteins"
VAR_DIR = "/data/projet1/02_Pipeline_Output/3_Mutations"
OUT_DIR = "/data/projet1/02_Pipeline_Output/8_Stats_Mutations"

# ---------- 1. proteines de souche avec codon stop premature ----------
ref_len = {}
cur = None
with open(REF) as fh:
    for l in fh:
        l = l.strip()
        if l.startswith(">"):
            cur = l[1:].split()[0]
            ref_len[cur] = 0
        elif cur:
            ref_len[cur] += len(l)

premature = {}   # fichier -> (nb, nb_sequences)
for pf in sorted(glob.glob(os.path.join(PROT_DIR, "*.prot.fasta"))):
    base = os.path.basename(pf)
    n_seq = 0
    n_prem = 0
    with open(pf) as fh:
        for l in fh:
            l = l.strip()
            if l.startswith(">"):
                n_seq += 1
            elif l:
                if "*" in l.rstrip("*"):
                    n_prem += 1
    if n_prem:
        premature[base] = (n_prem, n_seq)

with open(os.path.join(OUT_DIR, "codons_stop_prematures.tsv"), "w", newline="") as fh:
    fh.write("fichier_prot\tn_proteines_stop_premature\tn_proteines_total\tpct\n")
    for b, (np_, nt) in sorted(premature.items(), key=lambda x: -x[1][0]):
        fh.write(f"{b}\t{np_}\t{nt}\t{100*np_/nt:.2f}\n")

tot_prem = sum(v[0] for v in premature.values())
tot_seq = sum(v[1] for v in premature.values())
print(f"=== CODONS STOP PREMATURES ===")
print(f"Proteines avec stop premature : {tot_prem} / {tot_seq} ({100*tot_prem/tot_seq:.2f} %)")
print(f"Fichiers concernes : {len(premature)}")

# ---------- 2. quels CDS de souche n'ont AUCUN variant (prot equente) ----------
# on recupere les ids qui ont au moins une ligne dans un .var
var_ids = set()
for vf in glob.glob(os.path.join(VAR_DIR, "*.var")):
    with open(vf) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                var_ids.add(p[0])

# id de chaque sequence dans les .prot.fasta
no_variant = []
for pf in sorted(glob.glob(os.path.join(PROT_DIR, "*.prot.fasta")))[:6]:
    base = os.path.basename(pf)
    with open(pf) as fh:
        for l in fh:
            if l.startswith(">"):
                sid = l[1:].strip()
                if sid not in var_ids:
                    no_variant.append((base, sid))

print(f"\n=== CDS SANS AUCUN VARIANT ===")
print(f"Ids presents dans les .var : {len(var_ids)}")
print(f"Exemple de sequences sans variant (sur 6 fichiers testes) : {len(no_variant)}")
for b, s in no_variant[:8]:
    print(f"  {b}  {s[:70]}")

print("\nSorties dans", OUT_DIR)
