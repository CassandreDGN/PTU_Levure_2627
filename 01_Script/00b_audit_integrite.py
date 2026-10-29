#!/usr/bin/env python3
"""Audit cible : prouver qu'aucun fichier du groupe n'a ete modifie par moi.
Verification croisee : auteur (sonmez/degen/meguebli) + date de modification."""
import subprocess, os

def run(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return r.stdout.strip()

print("=== A. Fichiers du groupe et leurs dates (doivent etre <= 21 oct) ===")
print(run("ls -la --time-style=long-iso /data/projet1/01_Script/ /data/projet1/05_Documentation/ /data/projet1/README.md /data/projet1/cahier_de_labo.md 2>&1"))

print("\n=== B. Fichiers modifies le 25 oct (hors conda et hors les miens) ===")
print(run("find /data/projet1 -newermt '2026-10-25 00:00' -type f "
          "-not -path '*/.git/*' -not -path '*/06_Conda/*' -not -path '*/03_Archive/*' "
          "-not -path '*/conda/*' -not -user sonmez "
          "-printf '%TY-%Tm-%Td %TH:%TM  %u  %p\\n' 2>/dev/null | sort"))

print("\n=== C. Fichiers possedes par sonmez, hors conda (ce sont les miens) ===")
print(run("find /data/projet1 -user sonmez -type f "
          "-not -path '*/.git/*' -not -path '*/06_Conda/*' -not -path '*/03_Archive/*' "
          "-not -path '*/conda/*' "
          "-printf '%TY-%Tm-%Td %TH:%TM  %10s  %p\\n' 2>/dev/null | sort"))

print("\n=== D. Mes scripts jamais numerotes 06/07 (trou de numerotation) ===")
print(run("ls /data/projet1/01_Script/"))

print("\n=== E. Le dossier var_normalises a-t-il ete cree ? ===")
print(run("ls -la /data/projet1/02_Pipeline_Output/8_Stats_Mutations/var_normalises 2>&1 | head -3"))

print("\n=== F. Place disque ===")
print(run("df -h / /data | grep -v Filesystem"))
