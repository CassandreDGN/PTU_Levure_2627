#!/usr/bin/env python3
"""Faits bruts pour la presentation : verification des choix techniques reels."""
import os

def t(p, n=2000):
    return open(p, encoding="utf-8", errors="replace").read()[:n]

print("=== 02.1_blastP.py : colonnes de sortie ===")
print(t("/data/projet1/01_Script/02.1_blastP.py", 700))
print("\n=== 02.2_diamond.py : colonnes de sortie ===")
print(t("/data/projet1/01_Script/02.2_diamond.py", 900))
print("\n=== resume_blast.py (Beray, 11 oct) : que fait-il ? ===")
print(t("/data/projet1/04_Resultats_BlastN/resume_blast.py", 1963))
print("\n=== extract_mutations.py (Beray, 19 oct) : logique ===")
print(t("/data/projet1/04_Resultats_BlastN/extract_mutations.py", 2422))
print("\n=== check_stops.py (Beray, 19 oct) ===")
print(t("/data/projet1/04_Resultats_BlastN/check_stops.py", 742))
print("\n=== detect_troncatures.py (Beray, 19 oct) : entete seulement ===")
print(t("/data/projet1/04_Resultats_BlastN/detect_troncatures.py", 900))
