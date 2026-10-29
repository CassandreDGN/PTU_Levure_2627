#!/usr/bin/env python3
"""Qualite des variants (sanity check) :
 - substitutions vers '*' (STOP)
 - deletions/insertions sur '*'
 - distribution de taille des del/ins
 - variants dont la position depasse la longueur de la proteine de reference
"""
import os, re, glob
from collections import Counter

VAR_DIR = "/data/projet1/02_Pipeline_Output/3_Mutations"
REF = "/data/projet1/00_Data/ref_seq_proteique_fasta/orf_trans_all.fasta"
OUT = "/data/projet1/02_Pipeline_Output/8_Stats_Mutations/sanity_variants.tsv"

# longueurs de reference
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

RE_SUB = re.compile(r"^([A-Z*])(\d+)([A-Z*])$")
RE_DEL = re.compile(r"^([A-Z*])(\d+)del$")
RE_INS = re.compile(r"^(\d+)_(\d+)ins([A-Z*]+)$")

sub_stop = 0          # substitutions dont l'AA mutant est un STOP
del_stop = 0          # deletions d'un STOP (douteux)
ins_stop = 0          # insertions contenant un STOP
sub_total = del_total = ins_total = 0
pos_out = 0           # position > longueur de la proteine de reference
ins_len = Counter()
del_suspicious = 0

files = sorted(glob.glob(os.path.join(VAR_DIR, "*.var")))
for vf in files:
    with open(vf) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            gene, var = p[0], p[1]
            m = RE_SUB.match(var)
            if m:
                sub_total += 1
                if m.group(3) == "*":
                    sub_stop += 1
                pos = int(m.group(2))
                if ref_len.get(gene, 0) and pos > ref_len[gene]:
                    pos_out += 1
                continue
            m = RE_DEL.match(var)
            if m:
                del_total += 1
                if m.group(1) == "*":
                    del_stop += 1
                pos = int(m.group(2))
                if ref_len.get(gene, 0) and pos > ref_len[gene]:
                    pos_out += 1
                continue
            m = RE_INS.match(var)
            if m:
                ins_total += 1
                seq = m.group(3)
                ins_len[len(seq)] += 1
                if "*" in seq:
                    ins_stop += 1
                pos = int(m.group(2))
                if ref_len.get(gene, 0) and pos > ref_len[gene]:
                    pos_out += 1

with open(OUT, "w") as fh:
    fh.write("metrique\tvaleur\n")
    fh.write(f"substitutions_total\t{sub_total}\n")
    fh.write(f"substitutions_vers_STOP\t{sub_stop}\n")
    fh.write(f"pct_sub_vers_STOP\t{100*sub_stop/max(sub_total,1):.2f}\n")
    fh.write(f"deletions_total\t{del_total}\n")
    fh.write(f"deletions_de_STOP\t{del_stop}\n")
    fh.write(f"insertions_total\t{ins_total}\n")
    fh.write(f"insertions_contenant_STOP\t{ins_stop}\n")
    fh.write(f"positions_hors_proteine_ref\t{pos_out}\n")

print("=== SANITY CHECK ===")
print(f"substitutions totales          : {sub_total}")
print(f"  dont vers STOP (*)           : {sub_stop}  ({100*sub_stop/max(sub_total,1):.2f} %)")
print(f"deletions totales              : {del_total}")
print(f"  dont deletion d'un STOP      : {del_stop}")
print(f"insertions totales             : {ins_total}")
print(f"  dont insertion avec STOP     : {ins_stop}")
print(f"positions > longueur ref       : {pos_out}")
print("\n=== distribution de taille des insertions (AA) ===")
for k in sorted(ins_len)[:10]:
    print(f"  {k} AA : {ins_len[k]}")
print(f"  ... max = {max(ins_len) if ins_len else 0} AA")
print("\nSortie :", OUT)
