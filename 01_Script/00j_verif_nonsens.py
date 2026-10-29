#!/usr/bin/env python3
"""Peut-on vraiment faire sens / non-sens ? Verification sur les CDS bruts.
Un codon non-sens = codon stop (TAA/TAG/TGA) en pleine sequence."""
import os, glob

CODONS_STOP = {"TAA", "TAG", "TGA"}

def lire_fasta(chemin, limite=None):
    """Renvoie {id: sequence}."""
    seqs, nom, buf = {}, None, []
    with open(chemin) as fh:
        for ligne in fh:
            ligne = ligne.strip()
            if ligne.startswith(">"):
                if nom is not None:
                    seqs[nom] = "".join(buf)
                    if limite and len(seqs) >= limite:
                        return seqs
                nom = ligne[1:].split()[0]
                buf = []
            elif nom is not None:
                buf.append(ligne)
    if nom is not None:
        seqs[nom] = "".join(buf)
    return seqs

print("=== 1. Le script 01_CDStoPROTEIN.py filtre-t-il les STOP prematures ? ===")
src = open("/data/projet1/01_Script/01_CDStoPROTEIN.py").read()
for mot in ["premature", "STOP", "stop", "rejec", "reject"]:
    nb = src.count(mot)
    if nb:
        print(f"  le mot '{mot}' apparait {nb} fois")
# montrer les lignes concernées
for i, l in enumerate(src.splitlines()):
    if "premature" in l.lower() or ("stop" in l.lower() and ("*" in l or "rejec" in l.lower())):
        print(f"    ligne {i+1}: {l.strip()[:110]}")

print("\n=== 2. Les CDS de souche contiennent-ils des codons stop internes ? ===")
# on prend quelques fichiers representatifs
fichiers = sorted(glob.glob("/data/projet1/00_Data/cds_fasta/*.cds.fa"))[:6]
print(f"(test sur {len(fichiers)} fichiers)\n")
for f in fichiers:
    seqs = lire_fasta(f)
    n_seq = len(seqs)
    n_stop = 0
    exemples = []
    for nom, seq in seqs.items():
        codons = [seq[i:i+3] for i in range(0, len(seq) - 2, 3)]
        # on ignore le dernier codon (stop terminal normal)
        internes = [c for c in codons[:-1] if c in CODONS_STOP]
        if internes:
            n_stop += 1
            if len(exemples) < 2:
                exemples.append((nom, len(seq), len(internes), internes[:3]))
    print(f"  {os.path.basename(f)[:45]:47s} {n_seq:5d} CDS, dont {n_stop:4d} avec stop interne ({100*n_stop/max(n_seq,1):.2f} %)")
    for nom, lg, nb, cod in exemples:
        print(f"      ex: {nom[:40]:42s} {lg} nt, {nb} stop interne(s) {cod}")

print("\n=== 3. Reference S288C : y a-t-il des stop internes (controle) ? ===")
ref = lire_fasta("/data/projet1/00_Data/reference_S288C/orf_coding_all.fasta")
n_stop_ref = 0
for nom, seq in ref.items():
    codons = [seq[i:i+3] for i in range(0, len(seq) - 2, 3)]
    if any(c in CODONS_STOP for c in codons[:-1]):
        n_stop_ref += 1
print(f"  {len(ref)} genes, dont {n_stop_ref} avec stop interne")
