#!/usr/bin/env python3
"""Listes de controle pre/post-intervention : hash + date de TOUS les fichiers du
groupe (hors donnees lourdes et conda). C'est la liste de reference pour verifier
que rien n'a change."""
import os, hashlib, subprocess

SKIP_DIRS = ("/.git/", "/06_Conda/", "/03_Archive/", "/conda/", "/00_Data/cds_fasta/")
ROOT = "/data/projet1"

out = []
for dirpath, dirnames, filenames in os.walk(ROOT):
    if any(s in dirpath for s in SKIP_DIRS):
        dirnames[:] = []
        continue
    for fn in sorted(filenames):
        p = os.path.join(dirpath, fn)
        try:
            st = os.stat(p)
        except OSError:
            continue
        if st.st_size > 200_000_000:
            out.append(f"{p}\t{st.st_size}\t{st.st_mtime_ns}\tBIG")
            continue
        try:
            h = hashlib.md5(open(p, "rb").read()).hexdigest()
        except OSError:
            h = "NA"
        out.append(f"{p}\t{st.st_size}\t{st.st_mtime_ns}\t{h}")

dest = "/data/projet1/02_Pipeline_Output/8_Stats_Mutations/manifeste_fichiers.txt"
with open(dest, "w") as fh:
    fh.write("\n".join(out) + "\n")
print(f"{len(out)} fiches -> {dest}")
print("Total fichiers hors donnees lourdes :", len(out))
