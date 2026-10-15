#!/usr/bin/env python3

import re
from pathlib import Path
from collections import defaultdict

# --- CONFIGURATION DES CHEMINS ---
VAR_DIR = Path("/data/projet1/02_Pipeline_Output/3_Mutations")
REF_FASTA = Path("/data/projet1/00_Data/ref_seq_proteique_fasta/orf_trans_all.fasta")
BASE_OUT_DIR = Path("/data/projet1/02_Pipeline_Output/5_Blosum_Scores")
BASE_OUT_DIR.mkdir(parents=True, exist_ok=True)

# --- SEUILS DE FONCTIONNALITÉ ---
# Perte relative maximale tolérée (ex: 0.25 = 25% de perte)
RELATIVE_LOSS_THRESHOLD = 0.25 
GAP_OPEN = -3
GAP_EXTEND = -1

# --- MATRICE BLOSUM62 COMPLÈTE (20x20) ---
BLOSUM62 = {
    'A': {'A': 4, 'R': -1, 'N': -2, 'D': -2, 'C': 0, 'Q': -1, 'E': -1, 'G': 0, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 0, 'W': -3, 'Y': -2, 'V': 0},
    'R': {'A': -1, 'R': 5, 'N': 0, 'D': -2, 'C': -3, 'Q': 1, 'E': 0, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 2, 'M': -1, 'F': -3, 'P': -2, 'S': -1, 'T': -1, 'W': -3, 'Y': -2, 'V': -3},
    'N': {'A': -2, 'R': 0, 'N': 6, 'D': 1, 'C': -3, 'Q': 0, 'E': 0, 'G': 0, 'H': 1, 'I': -3, 'L': -3, 'K': 0, 'M': -2, 'F': -3, 'P': -2, 'S': 1, 'T': 0, 'W': -4, 'Y': -2, 'V': -3},
    'D': {'A': -2, 'R': -2, 'N': 1, 'D': 6, 'C': -3, 'Q': 0, 'E': 2, 'G': -1, 'H': -1, 'I': -3, 'L': -4, 'K': -1, 'M': -3, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3},
    'C': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': 9, 'Q': -3, 'E': -4, 'G': -3, 'H': -3, 'I': -1, 'L': -1, 'K': -3, 'M': -1, 'F': -2, 'P': -3, 'S': -1, 'T': -1, 'W': -2, 'Y': -2, 'V': -1},
    'Q': {'A': -1, 'R': 1, 'N': 0, 'D': 0, 'C': -3, 'Q': 5, 'E': 2, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 1, 'M': 0, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -2, 'Y': -1, 'V': -2},
    'E': {'A': -1, 'R': 0, 'N': 0, 'D': 2, 'C': -4, 'Q': 2, 'E': 5, 'G': -2, 'H': 0, 'I': -3, 'L': -3, 'K': 1, 'M': -2, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2},
    'G': {'A': 0, 'R': -2, 'N': 0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G': 6, 'H': -2, 'I': -4, 'L': -4, 'K': -2, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -2, 'W': -2, 'Y': -3, 'V': -3},
    'H': {'A': -2, 'R': 0, 'N': 1, 'D': -1, 'C': -3, 'Q': 0, 'E': 0, 'G': -2, 'H': 8, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -1, 'P': -2, 'S': -1, 'T': -2, 'W': -2, 'Y': 2, 'V': -3},
    'I': {'A': -1, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -3, 'E': -3, 'G': -4, 'H': -3, 'I': 4, 'L': 2, 'K': -3, 'M': 1, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -3, 'Y': -1, 'V': 3},
    'L': {'A': -1, 'R': -2, 'N': -3, 'D': -4, 'C': -1, 'Q': -2, 'E': -3, 'G': -4, 'H': -3, 'I': 2, 'L': 4, 'K': -2, 'M': 2, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -1, 'V': 1},
    'K': {'A': -1, 'R': 2, 'N': 0, 'D': -1, 'C': -3, 'Q': 1, 'E': 1, 'G': -2, 'H': -1, 'I': -3, 'L': -2, 'K': 5, 'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2},
    'M': {'A': -1, 'R': -1, 'N': -2, 'D': -3, 'C': -1, 'Q': 0, 'E': -2, 'G': -3, 'H': -2, 'I': 1, 'L': 2, 'K': -1, 'M': 5, 'F': 0, 'P': -2, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': 1},
    'F': {'A': -2, 'R': -3, 'N': -3, 'D': -3, 'C': -2, 'Q': -3, 'E': -3, 'G': -3, 'H': -1, 'I': 0, 'L': 0, 'K': -3, 'M': 0, 'F': 6, 'P': -4, 'S': -2, 'T': -2, 'W': 1, 'Y': 3, 'V': -1},
    'P': {'A': -1, 'R': -2, 'N': -2, 'D': -1, 'C': -3, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -4, 'P': 7, 'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -2},
    'S': {'A': 1, 'R': -1, 'N': 1, 'D': 0, 'C': -1, 'Q': 0, 'E': 0, 'G': 0, 'H': -1, 'I': -2, 'L': -2, 'K': 0, 'M': -1, 'F': -2, 'P': -1, 'S': 4, 'T': 1, 'W': -3, 'Y': -2, 'V': -2},
    'T': {'A': 0, 'R': -1, 'N': 0, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 5, 'W': -2, 'Y': -2, 'V': 0},
    'W': {'A': -3, 'R': -3, 'N': -4, 'D': -4, 'C': -2, 'Q': -2, 'E': -3, 'G': -2, 'H': -2, 'I': -3, 'L': -2, 'K': -3, 'M': -1, 'F': 1, 'P': -4, 'S': -3, 'T': -2, 'W': 11, 'Y': 2, 'V': -3},
    'Y': {'A': -2, 'R': -2, 'N': -2, 'D': -3, 'C': -2, 'Q': -1, 'E': -2, 'G': -3, 'H': 2, 'I': -1, 'L': -1, 'K': -2, 'M': -1, 'F': 3, 'P': -3, 'S': -2, 'T': -2, 'W': 2, 'Y': 7, 'V': -1},
    'V': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -2, 'E': -2, 'G': -3, 'H': -3, 'I': 3, 'L': 1, 'K': -2, 'M': 1, 'F': -1, 'P': -2, 'S': -2, 'T': 0, 'W': -3, 'Y': -1, 'V': 4}
}

def get_blosum_score(aa1, aa2):
    aa1_u, aa2_u = aa1.upper(), aa2.upper()
    if aa1_u in BLOSUM62 and aa2_u in BLOSUM62[aa1_u]:
        return BLOSUM62[aa1_u][aa2_u]
    return -4  # Pénalité par défaut si acide aminé inconnu (ex: X)

def parse_variant(var_str):
    var_str = var_str.strip()
    sub_match = re.match(r'^([A-Z])([0-9]+)([A-Z])$', var_str)
    if sub_match:
        return 'substitution', int(sub_match.group(2)) - 1, sub_match.group(1), sub_match.group(3)
    if 'del' in var_str:
        numbers = [int(n) for n in re.findall(r'[0-9]+', var_str)]
        pos = max(0, numbers[0] - 1) if numbers else 0
        length = (abs(numbers[1] - numbers[0]) + 1) if len(numbers) >= 2 else 1
        return 'deletion', pos, length, None
    if 'ins' in var_str:
        numbers = [int(n) for n in re.findall(r'[0-9]+', var_str)]
        pos = max(0, numbers[0] - 1) if numbers else 0
        ins_seq = var_str.split('ins')[-1]
        length = len(ins_seq) if ins_seq.isalpha() else 1
        return 'insertion', pos, length, ins_seq
    return 'unknown', 0, 0, None

# 1. Chargement du FASTA et calcul des Scores Sauvages
print("Chargement du FASTA et calcul des scores sauvages de référence...")
proteins = {}
wt_scores = {}
current_id = None
current_seq = []

with open(REF_FASTA) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current_id:
                seq = "".join(current_seq)
                proteins[current_id] = seq
                wt_scores[current_id] = sum(get_blosum_score(aa, aa) for aa in seq)
            current_id = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line)
    if current_id:
        seq = "".join(current_seq)
        proteins[current_id] = seq
        wt_scores[current_id] = sum(get_blosum_score(aa, aa) for aa in seq)

print(f"-> {len(proteins)} protéines de référence chargées avec leurs scores parfaits.")

var_files = list(VAR_DIR.glob("*.var"))
print(f"Traitement de {len(var_files)} fichiers de variants...")

for var_file in sorted(var_files):
    copy_name = var_file.name.replace('.var', '')
    strain_dir = BASE_OUT_DIR / copy_name
    strain_dir.mkdir(parents=True, exist_ok=True)
    
    prot_variants = defaultdict(list)
    with open(var_file) as fin:
        for line in fin:
            parts = line.strip().split()
            if len(parts) >= 2:
                prot_variants[parts[0]].append(parts[1])
                
    gene_loss = defaultdict(int)
    gene_mut_counts = defaultdict(int)
    gene_mut_list = defaultdict(list)
    
    with open(strain_dir / "variant_details.tsv", "w") as fout:
        fout.write("protein\tvariant\ttype\tscore_delta\tstatus\n")
        
        for prot_id, var_strs in prot_variants.items():
            if prot_id not in proteins:
                continue
            seq = proteins[prot_id]
            
            parsed = []
            for v in var_strs:
                v_type, pos, ref_or_len, alt = parse_variant(v)
                parsed.append((pos, v_type, ref_or_len, alt, v))
                
            parsed.sort(key=lambda x: x[0])
            
            grouped = []
            curr = []
            for p in parsed:
                if not curr:
                    curr = [p]
                else:
                    last = curr[-1]
                    if p[1] == last[1] == 'deletion' and p[0] == last[0] + last[2]:
                        curr.append(p)
                    elif p[1] == last[1] == 'insertion' and p[0] == last[0]:
                        curr.append(p)
                    else:
                        grouped.append(curr)
                        curr = [p]
            if curr:
                grouped.append(curr)

            for grp in grouped:
                v_type = grp[0][1]
                raw_strs = ",".join([x[4] for x in grp])
                delta_score = 0
                status = "OK"
                
                if v_type == 'substitution':
                    pos, ref, alt = grp[0][0], grp[0][2], grp[0][3]
                    if 0 <= pos < len(seq):
                        expected_ref = seq[pos]
                        # La perte = (Score de la mutation) - (Score parfait de l'acide aminé d'origine)
                        delta_score = get_blosum_score(expected_ref, alt) - get_blosum_score(expected_ref, expected_ref)
                        if expected_ref != ref:
                            status = f"REF_MISMATCH"
                    else:
                        delta_score = -4
                        status = "OUT_OF_BOUNDS"
                elif v_type in ('deletion', 'insertion'):
                    total_length = sum([x[2] for x in grp])
                    delta_score = GAP_OPEN + (total_length - 1) * GAP_EXTEND
                else:
                    delta_score = -2
                    status = "UNPARSED"
                
                fout.write(f"{prot_id}\t{raw_strs}\t{v_type}\t{delta_score}\t{status}\n")
                
                gene_loss[prot_id] += delta_score  # S'additionne en négatif
                gene_mut_counts[prot_id] += 1
                gene_mut_list[prot_id].append(raw_strs)

    # Synthèse avec la formule de perte relative
    with open(strain_dir / "gene_summary.tsv", "w") as fgene:
        fgene.write("protein\tscore_sauvage\tscore_mute\tdelta_relatif\tn_mutations\tmutations_list\tis_functional\n")
        
        for prot_id in sorted(proteins.keys()):
            wt_score = wt_scores[prot_id]
            loss = gene_loss[prot_id] # Valeur négative
            mutated_score = wt_score + loss
            
            # Calcul du delta relatif (perte) = (WT - Mute) / WT 
            delta_relatif = (wt_score - mutated_score) / wt_score if wt_score > 0 else 0
            
            n_mut = gene_mut_counts[prot_id]
            muts = ";".join(gene_mut_list[prot_id]) if n_mut > 0 else "WT"
            
            # Décision finale : Inactif si la perte relative dépasse le seuil (ex: > 0.25)
            if delta_relatif > RELATIVE_LOSS_THRESHOLD:
                is_func = "NO"
            else:
                is_func = "YES"
                
            fgene.write(f"{prot_id}\t{wt_score}\t{mutated_score}\t{delta_relatif:.4f}\t{n_mut}\t{muts}\t{is_func}\n")

print("\nÉtape 3 ou 4 terminée : Synthèses de perte relative générées avec succès.")