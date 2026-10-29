#!/usr/bin/env python3
"""Le pipeline actuel decide-il jamais NON fonctionnel ? Et que contiendrait
vraiment une classification sens/faux-sens ?"""
import glob, os

print("=== 1. Toutes les valeurs possibles de is_functional ===")
r = os.popen("cat /data/projet1/02_Pipeline_Output/5_Blosum_Scores/*/gene_summary.tsv 2>/dev/null | cut -f7 | sort | uniq -c | sort -rn").read()
print(r)

print("=== 2. delta_relatif : max et distribution ===")
print("max :")
print(os.popen("cat /data/projet1/02_Pipeline_Output/5_Blosum_Scores/*/gene_summary.tsv 2>/dev/null | awk -F'\\t' 'NR>1 && $4!=\"\" {print $4}' | sort -g | tail -5").read())
print("combien de genes ont delta_relatif > 0.25 (seuil declare dans 03.2) :")
print(os.popen("cat /data/projet1/02_Pipeline_Output/5_Blosum_Scores/*/gene_summary.tsv 2>/dev/null | awk -F'\\t' 'NR>1 && $4+0 > 0.25' | wc -l").read())
print("combien de genes ont n_mutations > 0 :")
print(os.popen("cat /data/projet1/02_Pipeline_Output/5_Blosum_Scores/*/gene_summary.tsv 2>/dev/null | awk -F'\\t' 'NR>1 && $5+0 > 0' | wc -l").read())
print("total de lignes :")
print(os.popen("cat /data/projet1/02_Pipeline_Output/5_Blosum_Scores/*/gene_summary.tsv 2>/dev/null | wc -l").read())

print("=== 3. variant_details.tsv : format reel ===")
f = sorted(glob.glob("/data/projet1/02_Pipeline_Output/5_Blosum_Scores/*/variant_details.tsv"))
if f:
    print("fichier :", os.path.basename(os.path.dirname(f[0])))
    print(open(f[0]).read()[:700])

print("=== 4. Le CDS de reference S288C contient-il les positions ? (utile pour le codon) ===")
print(os.popen("head -1 /data/projet1/00_Data/reference_S288C/orf_coding_all.fasta").read())
print(os.popen("grep -c '^>' /data/projet1/00_Data/reference_S288C/orf_coding_all.fasta").read())

print("=== 5. Pour YAL001C, ou tombe la position 432 dans le codon ? ===")
# charger la sequence d'YAL001C et examiner les codons autour de la position 432
seq = ""
name = None
for ligne in open("/data/projet1/00_Data/reference_S288C/orf_coding_all.fasta"):
    ligne = ligne.strip()
    if ligne.startswith(">"):
        if name == "YAL001C":
            break
        name = ligne[1:].split()[0]
        seq = ""
    else:
        seq += ligne
if name == "YAL001C":
    print("longueur YAL001C (nt) :", len(seq), "->", len(seq)//3, "codons")
    for pos_prot in [432, 435, 567, 731, 1089]:
        i = (pos_prot - 1) * 3
        codon = seq[i:i+3]
        cadre = "en phase" if i % 3 == 0 else "HORS PHASE"
        print(f"  position proteine {pos_prot} -> nucleotides {i+1}-{i+3} : codon {codon} ({cadre})")
else:
    print("YAL001C non trouve (dernier nom lu :", name, ")")
