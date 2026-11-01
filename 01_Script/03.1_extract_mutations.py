#!/usr/bin/env python3
"""
Étape 3 : extraction des variants 
Modifié pour SIFT4G : 
- Extrait uniquement les substitutions pour les fichiers .subst.
- Génère un fichier `loss_of_function_events.tsv` (STOP prématurés et Délétions >= 30% de la prot).
"""

import os
import glob
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import Align

STAT_COLS = [
    "n_hits", "n_ref_genes", "collisions", "dropped_collisions",
    "n_mutations", "sub", "skipped_del", "skipped_ins", "skipped_stop",
    "realigned", "skipped_X", "skipped_length_mismatch", "skipped_malformed"
] 

def extract_mutations_subs_only(sseq, qseq, sseqid, lof_file, ref_start=1):
    """
    Parcourt un alignement gappé et extrait UNIQUEMENT les substitutions.
    Évalue les STOP prématurés et les grandes délétions (>= 30%) pour le fichier LoF.
    """
    subs = []
    stats = Counter()
    ref_pos = ref_start - 1
    
    # Calcul de la longueur de la protéine de référence (sans les gaps d'alignement)
    ref_len = len(sseq.replace("-", ""))
    total_deleted_aa = 0

    for r, q in zip(sseq, qseq):
        # Mise à jour de la position sur la référence
        if r != '-':
            ref_pos += 1

        # Ignorer les X
        if r == 'X' or q == 'X':
            stats["skipped_X"] += 1
            continue

        # Détection des STOP prématurés (gain de STOP dans la requête)
        if q == '*' and r != '*':
            stats["skipped_stop"] += 1
            lof_file.write(f"{sseqid}\tPREMATURE_STOP\t{r}{ref_pos}{q}\n")
            continue
        elif r == '*':
            stats["skipped_stop"] += 1
            continue

        # Détection des variants
        if r != q:
            if r == '-':
                stats["skipped_ins"] += 1
                # Les insertions ne sont plus reportées comme LoF
            elif q == '-':
                stats["skipped_del"] += 1
                total_deleted_aa += 1
            else:
                # C'est une substitution pure (ex: R42K)
                subs.append(f"{r}{ref_pos}{q}")
                stats["sub"] += 1

    # Évaluation de la délétion globale à la fin de la séquence
    if ref_len > 0 and (total_deleted_aa / ref_len) >= 0.30:
        lof_file.write(f"{sseqid}\tLARGE_DELETION_>30%\t{total_deleted_aa}_aa_lost\n")

    return subs, stats

def process_tsv_single(args_tuple):
    tsv_file, out_dir, best_per_gene = args_tuple
    
    filename = os.path.basename(tsv_file)
    strain_name = filename.replace(".fasta_blast.tsv", "").replace("_blast.tsv", "")
    
    strain_dir = os.path.join(out_dir, strain_name)
    subst_dir = os.path.join(strain_dir, "subst")
    os.makedirs(subst_dir, exist_ok=True)
    
    lof_path = os.path.join(strain_dir, "loss_of_function_events.tsv")
    
    stats = Counter()
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    try:
        hits = []
        with open(tsv_file, "r") as f_tsv:
            for line in f_tsv:
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) < 8:
                    stats["skipped_malformed"] += 1
                    continue
                hits.append(parts)

        stats["n_hits"] = len(hits)
        stats["n_ref_genes"] = len({p[1] for p in hits})
        stats["collisions"] = stats["n_hits"] - stats["n_ref_genes"]

        if best_per_gene:
            best = {}
            for p in hits:
                if p[1] not in best or float(p[5]) > float(best[p[1]][5]):
                    best[p[1]] = p
            stats["dropped_collisions"] = len(hits) - len(best)
            hits = list(best.values())

        gene_to_subs = {}
        
        with open(lof_path, "w") as f_lof:
            f_lof.write("protein\tevent_type\tposition\n")
            
            for p in hits:
                sseqid, qseq, sseq = p[1], p[6], p[7]
                try:
                    sstart = int(p[8]) if len(p) > 8 else 1
                except ValueError:
                    stats["skipped_malformed"] += 1
                    continue

                if qseq == sseq:
                    continue                                   

                if len(sseq) != len(qseq) and '-' not in sseq and '-' not in qseq:
                    aln = aligner.align(sseq, qseq)[0]
                    sseq, qseq = str(aln[0]), str(aln[1])
                    stats["realigned"] += 1

                if len(sseq) != len(qseq):                    
                    stats["skipped_length_mismatch"] += 1
                    continue

                subs, st = extract_mutations_subs_only(sseq, qseq, sseqid, f_lof, ref_start=sstart)
                stats.update(st)

                if subs:
                    if sseqid not in gene_to_subs:
                        gene_to_subs[sseqid] = []
                    gene_to_subs[sseqid].extend(subs)

        stats["n_mutations"] = sum(len(v) for v in gene_to_subs.values())

        for gene_id, sub_list in gene_to_subs.items():
            unique_subs = list(set(sub_list)) 
            with open(os.path.join(subst_dir, f"{gene_id}.subst"), "w") as f_sub:
                f_sub.write("\n".join(unique_subs) + "\n")

        return strain_name, dict(stats), None
    except Exception as e:
        return strain_name, dict(stats), str(e)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extracts ONLY substitutions and formats for SIFT4G. Logs LoF events (Stops & >30% Deletions).")
    parser.add_argument("-t", "--tsv_dir", default="/data/projet1/02_Pipeline_Output/2_Blast")
    parser.add_argument("-o", "--out_dir", default="/data/projet1/02_Pipeline_Output/3_SIFT_Mutations")
    parser.add_argument("-j", "--jobs", type=int, default=8)
    parser.add_argument("--best-per-gene", action="store_true")

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    tsv_files = sorted(glob.glob(os.path.join(args.tsv_dir, "*.tsv")))
    if not tsv_files:
        print(f"Error: No TSV files found in {args.tsv_dir}", flush=True)
        exit(1)

    print(f"Running extraction on {len(tsv_files)} files with {args.jobs} cores...", flush=True)

    tasks = [(tsv, args.out_dir, args.best_per_gene) for tsv in tsv_files]
    completed = 0
    errors = 0
    total = Counter()
    per_file = {}

    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(process_tsv_single, task) for task in tasks]
        for future in as_completed(futures):
            strain_name, st, err = future.result()
            if err:
                errors += 1
                print(f"[{completed + errors}/{len(tsv_files)}] {strain_name} -> ERROR: {err}", flush=True)
            else:
                completed += 1
                total.update(st)
                per_file[strain_name] = st

            if (completed + errors) % 50 == 0 or (completed + errors) == len(tsv_files):
                print(f"Progress: {completed + errors}/{len(tsv_files)} files processed...", flush=True)

    stats_path = os.path.join(args.out_dir, "_extraction_stats.tsv")
    with open(stats_path, "w") as f_stats:
        f_stats.write("file\t" + "\t".join(STAT_COLS) + "\n")
        for name in sorted(per_file):
            f_stats.write(name + "\t" + "\t".join(str(per_file[name].get(c, 0)) for c in STAT_COLS) + "\n")

    print(f"\nStep 3 completed: Structure generated in {args.out_dir} "
          f"({total['sub']} substitutions extracted. {total['skipped_del']+total['skipped_ins']} indels and {total['skipped_stop']} stops recorded. {errors} errors).", flush=True)