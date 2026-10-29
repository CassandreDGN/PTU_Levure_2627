#!/usr/bin/env python3
"""Manifeste leger : hash des SEULS fichiers sources du groupe (scripts, docs, README).
Pas de donnees lourdes. Ecrit le manifeste sur le serveur."""
import os, hashlib

ROOT = "/data/projet1"
CIBLES = [
    "01_Script",
    "05_Documentation",
    "README.md",
]
dest = os.path.join(ROOT, "02_Pipeline_Output/8_Stats_Mutations/manifeste_sources.txt")

out = []
for cible in CIBLES:
    chemin = os.path.join(ROOT, cible)
    fichiers = []
    if os.path.isfile(chemin):
        fichiers = [chemin]
    else:
        for dp, dn, fn in os.walk(chemin):
            for f in sorted(fn):
                fichiers.append(os.path.join(dp, f))
    for p in sorted(fichiers):
        try:
            data = open(p, "rb").read()
            h = hashlib.md5(data).hexdigest()
            st = os.stat(p)
            out.append(f"{h}  {st.st_size}  {p}")
        except OSError as e:
            out.append(f"ERREUR  {p}  {e}")

with open(dest, "w") as fh:
    fh.write("\n".join(out) + "\n")
print(f"{len(out)} fichiers sources hash es")
print("->", dest)
for l in out:
    print(l)
