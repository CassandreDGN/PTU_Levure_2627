#!/usr_bin/env python3
"""
1_translate_and_qc.py
---------------------

"""

import argparse
import os
import sys

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

VALID_NUCLEOTIDES = set("ACGT") #to make sure the sequence are made of valid nucleotides


def parse_fasta(fasta_path):
    #lit le fichier fasta et retourne un dictionnaire {header: sequence}
    sequences = {}
    current_header = None
    current_seq = []
    seen_headers = set()

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header:
                    seq_str = "".join(current_seq)
                    sequences[current_header] = seq_str
                current_header = line[1:].split()[0]  # Garde l'ID principal
                
                # Gestion des doublons
                if current_header in seen_headers:
                    print(f"  [WARNING] Header en doublon détecté : {current_header}")
                seen_headers.add(current_header)
                current_seq = []
            else:
                current_seq.append(line.upper())

        if current_header:
            sequences[current_header] = "".join(current_seq)

    return sequences


def analyze_and_translate_cds(cds_seq, min_length):
    """
    Analyse la qualité du CDS et le traduit si valide.
    
    Retourne :
        - status (str) : "VALID" ou le motif du rejet
        - result (str) : Séquence protéique si VALID, sinon détail/position du rejet
    """
    cds_seq = cds_seq.upper().strip()

    # Première verif : Séquence vide
    if not cds_seq:
        return "REJECTED_EMPTY_SEQUENCE", "Séquence nulle"

    # Deuxième verif : Longueur multiple de 3, si pas le cas on peut pas traduire toute façon
    if len(cds_seq) % 3 != 0:
        return "REJECTED_NOT_MULTIPLE_OF_3", f"Longueur = {len(cds_seq)} nt"

    # Troisième verif : Vérification du codon start ATG
    start_codon = cds_seq[:3]
    if start_codon != "ATG":
        return "REJECTED_NO_START_ATG", f"Codon initial : {start_codon}"

    # Quatrième verif : Traduction et vérification des codons
    protein = []
    codons = [cds_seq[i:i+3] for i in range(0, len(cds_seq), 3)]
    num_codons = len(codons)

    for idx, codon in enumerate(codons):
        # Vérification si le codon contient des nucléotides inconnus/invalides
        if not set(codon).issubset(VALID_NUCLEOTIDES) or codon not in CODON_TABLE:
            pos_nt = (idx * 3) + 1
            return "REJECTED_INVALID_CODON", f"Codon inconnu '{codon}' à la position nt {pos_nt}"

        aa = CODON_TABLE[codon]

        # Vérification des codons STOP prématurés
        if aa == "*":
            if idx < num_codons - 1:
                pos_nt = (idx * 3) + 1
                return "REJECTED_PREMATURE_STOP", f"Codon STOP '{codon}' à la position nt {pos_nt} (codon {idx+1}/{num_codons})"
            else:
                # Codon STOP normal en fin de séquence
                break

        protein.append(aa)

    protein_seq = "".join(protein)

    # Cinquième verif : Vérification de la longueur minimale de la protéine (paramètre ajustabe) 
    if len(protein_seq) < min_length:
        return "REJECTED_TOO_SHORT", f"Taille = {len(protein_seq)} AA (seuil : {min_length} AA)"

    return "VALID", protein_seq


def process_file(input_file, output_dir, min_length):
    """Traite un fichier FASTA unique."""
    base_name = os.path.basename(input_file)
    sample_id = os.path.splitext(base_name)[0]
    output_fasta_path = os.path.join(output_dir, f"{sample_id}.prot.fasta")

    sequences = parse_fasta(input_file)
    valid_proteins = {}
    rejected_records = []

    stats = {
        "sample": sample_id,
        "total": len(sequences),
        "valid": 0,
        "not_mult_3": 0,
        "no_atg": 0,
        "invalid_codon": 0,
        "premature_stop": 0,
        "too_short": 0,
        "empty": 0
    }

    for header, cds in sequences.items():
        status, detail = analyze_and_translate_cds(cds, min_length)

        if status == "VALID":
            valid_proteins[header] = detail
            stats["valid"] += 1
        else:
            rejected_records.append((sample_id, header, status, detail))
            if status == "REJECTED_NOT_MULTIPLE_OF_3":
                stats["not_mult_3"] += 1
            elif status == "REJECTED_NO_START_ATG":
                stats["no_atg"] += 1
            elif status == "REJECTED_INVALID_CODON":
                stats["invalid_codon"] += 1
            elif status == "REJECTED_PREMATURE_STOP":
                stats["premature_stop"] += 1
            elif status == "REJECTED_TOO_SHORT":
                stats["too_short"] += 1
            elif status == "REJECTED_EMPTY_SEQUENCE":
                stats["empty"] += 1

    # Écriture du FASTA protéique
    with open(output_fasta_path, "w") as out:
        for header, prot in valid_proteins.items():
            out.write(f">{header}\n{prot}\n")

    return stats, rejected_records


def main():
    parser = argparse.ArgumentParser(
        description="Étape 1 : Traduction CDS en protéines et Contrôle Qualité (QC) strict."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Chemin vers un fichier FASTA CDS ou un dossier contenant des fichiers FASTA"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Dossier de sortie pour les protéines et les rapports"
    )
    parser.add_argument(
        "-m", "--min-length", type=int, default=30,
        help="Seuil de longueur minimale en acides aminés (Défaut : 30)"
    )

    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    # Lister les fichiers d'entrée
    input_files = []
    if os.path.isfile(args.input):
        input_files.append(args.input)
    elif os.path.isdir(args.input):
        for root, _, files in os.walk(args.input):
            for file in sorted(files):
                if file.endswith((".fasta", ".fa", ".cds", ".txt")):
                    input_files.append(os.path.join(root, file))
    else:
        print(f"[ERROR] Entrée introuvable : {args.input}")
        sys.exit(1)

    print("=" * 60)
    print(" PIPELINE M2 - ÉTAPE 1 : TRADUCTION & CONTRÔLE QUALITÉ ")
    print("=" * 60)
    print(f"Fichiers à traiter : {len(input_files)}")
    print(f"Seuil de longueur  : {args.min_length} AA")
    print(f"Dossier de sortie  : {args.output}\n")

    all_stats = []
    all_rejected = []

    for idx, file_path in enumerate(input_files, 1):
        filename = os.path.basename(file_path)
        print(f"[{idx}/{len(input_files)}] Traitement de : {filename}")
        
        stats, rejected = process_file(file_path, args.output, args.min_length)
        all_stats.append(stats)
        all_rejected.extend(rejected)

        # Affichage interactif des résultats en console
        print(f"  ├─ CDS analysés        : {stats['total']}")
        print(f"  ├─ Protéines VALIDES   : {stats['valid']}")
        print(f"  ├─ Rejets (Non mult 3) : {stats['not_mult_3']}")
        print(f"  ├─ Rejets (Pas d'ATG)  : {stats['no_atg']}")
        print(f"  ├─ Rejets (Codon inv)  : {stats['invalid_codon']}")
        print(f"  ├─ Rejets (STOP préma) : {stats['premature_stop']}")
        print(f"  └─ Rejets (< {args.min_length} AA)    : {stats['too_short']}\n")

    # Écriture du rapport global de synthèse
    summary_path = os.path.join(args.output, "qc_summary.tsv")
    with open(summary_path, "w") as f:
        f.write("souche\ttotal_cds\tvalid_prot\tnot_multiple_3\tno_atg\tinvalid_codon\tpremature_stop\ttoo_short\tempty\n")
        for st in all_stats:
            f.write(
                f"{st['sample']}\t{st['total']}\t{st['valid']}\t{st['not_mult_3']}\t"
                f"{st['no_atg']}\t{st['invalid_codon']}\t{st['premature_stop']}\t"
                f"{st['too_short']}\t{st['empty']}\n"
            )

    # Écriture du rapport détaillé des séquences rejetées
    rejected_path = os.path.join(args.output, "rejected_sequences.tsv")
    with open(rejected_path, "w") as f:
        f.write("souche\tgene_id\tmotif_rejet\tdetail_position\n")
        for rej in all_rejected:
            f.write(f"{rej[0]}\t{rej[1]}\t{rej[2]}\t{rej[3]}\n")

    print("=" * 60)
    print("[DONE] TRAITEMENT TERMINÉ AVEC SUCCÈS")
    print(f"  ├─ Rapport de synthèse : {summary_path}")
    print(f"  └─ Rapport des rejets  : {rejected_path}")
    print("=" * 60)


if __name__ == "__main__":
    main()