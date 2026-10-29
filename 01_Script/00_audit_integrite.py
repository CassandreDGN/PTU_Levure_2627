#!/usr/bin/env python3
"""Audit : prouver qu'aucun fichier preexistant n'a ete modifie.
Liste ce qui est a moi (sonmez) et ce qui est intact, avec dates de modification."""
import subprocess, os

def run(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return r.stdout

print("=== 1. Fichiers appartenant a sonmez (les miens) ===")
print(run("find /data/projet1 -user sonmez -printf '%TY-%Tm-%Td %TH:%TM  %10s  %p\\n' 2>/dev/null | sort"))

print("=== 2. Dossier 8_Stats_Mutations (cree par moi) ===")
print(run("ls -la /data/projet1/02_Pipeline_Output/8_Stats_Mutations/"))

print("=== 3. Scripts que j'ai ecrits ===")
print(run("ls -la /data/projet1/01_Script/"))

print("=== 4. Derniere modif des fichiers du groupe (doit etre <= 25 oct) ===")
print(run("ls -la /data/projet1/01_Script/ /data/projet1/05_Documentation/ /data/projet1/README.md"))

print("=== 5. Rien de modifie apres le 25 oct 14:30 hors mes fichiers ===")
print(run("find /data/projet1 -newermt '2026-10-25 14:30' -not -user sonmez -not -path '*/.git/*' -printf '%TY-%Tm-%Td %TH:%TM  %u  %p\\n' 2>/dev/null | sort | head -40"))

print("=== 6. /data/projet1/tmp ===")
print(run("ls -la /data/projet1/tmp/ 2>&1 | head"))

print("=== 7. Etat git (rien de staging, rien de commite) ===")
print(run("cd /data/projet1 && git -c safe.directory=/data/projet1 status --short 2>&1 | head -20"))
print(run("cd /data/projet1 && git -c safe.directory=/data/projet1 log --oneline -3 2>&1"))
