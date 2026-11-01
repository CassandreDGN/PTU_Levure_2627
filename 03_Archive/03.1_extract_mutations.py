#!/usr/bin/env python3
"""
Étape 3 : extraction des variants à partir des alignements DIAMOND (outfmt 6).

Colonnes attendues dans les TSV (ordre de la commande de l'étape 2) :
  0 qseqid | 1 sseqid | 2 pident | 3 length | 4 evalue | 5 bitscore | 6 qseq | 7 sseq | 8 sstart

IMPORTANT : dans DIAMOND, les champs qseq/sseq sont la partie alignée SANS gaps.
Quand la souche et la référence ont des longueurs différentes (indel), l'alignement
est donc refait avec Biopython (global, gap open -10, gap extend -0.5) pour retrouver
les positions des gaps. Sans ce réalignement, tous les hits avec indel seraient perdus.

Sorties :
  <out_dir>/<souche>.var   : une ligne par variant -> sseqid<TAB>variant
                             (+ <TAB>qseqid si --with-qseqid)
                             AUCUN stop prématuré dedans (PROVEAN ne gère pas '*').
  <out_dir>/<souche>.stops : une ligne par stop prématuré -> sseqid<TAB>variant (ex. YAL001C  Y5*)
                             (+ <TAB>qseqid si --with-qseqid)
  <out_dir>/_extraction_stats.tsv : une ligne de statistiques par fichier

Notation (positions 1-based sur la référence S288C), inchangée par rapport à l'ancienne version :
  substitution        R42K
  délétion            R42del       (une ligne par résidu délété)
  insertion           42_43insV    (une ligne par résidu inséré)
  stop prématuré      Y5*          (fichier .stops uniquement) = présent dans un fichier à part
"""

import os
import glob
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

from Bio import Align

STAT_COLS = [
    "n_hits", "n_ref_genes", "collisions", "dropped_collisions",
    "n_mutations", "sub", "del", "ins", "stop_gain", "realigned",
    "skipped_X", "skipped_stop_loss", "skipped_stop_insertion",
    "skipped_length_mismatch", "skipped_malformed",
]


def extract_mutations(sseq, qseq, ref_start=1, truncate_at_stop=False):
    """
    Parcourt un alignement gappé (sseq = référence, qseq = souche, même longueur).
    Retourne (variants, stops prématurés, Counter de stats).

      - colonnes avec X : ignorées (comptées) ;
      - '*' dans qseq face à un résidu de la référence : stop prématuré, mis dans
        la liste des stops (jamais dans les variants) ;
      - '*' dans sseq (perte du stop de référence) : ignoré (compté).
    """
    muts, stops = [], []
    stats = Counter()
    ref_pos = ref_start - 1

    for r, q in zip(sseq, qseq):
        if r != '-':
            ref_pos += 1

        if r == 'X' or q == 'X':
            stats["skipped_X"] += 1
            continue

        if q == '*':
            if r == '-':
                stats["skipped_stop_insertion"] += 1
            else:
                if r != '*':                               # stop prématuré
                    stops.append(f"{r}{ref_pos}*")
                    stats["stop_gain"] += 1
                if truncate_at_stop:
                    break
            continue

        if r == '*':
            stats["skipped_stop_loss"] += 1
            continue

        if r != q:
            if r == '-':                                   # insertion
                muts.append(f"{ref_pos}_{ref_pos + 1}ins{q}")
                stats["ins"] += 1
            elif q == '-':                                 # délétion
                muts.append(f"{r}{ref_pos}del")
                stats["del"] += 1
            else:                                          # substitution
                muts.append(f"{r}{ref_pos}{q}")
                stats["sub"] += 1

    return muts, stops, stats


def process_tsv_single(args_tuple):
    tsv_file, out_dir, best_per_gene, with_qseqid, truncate_at_stop = args_tuple
    filename = os.path.basename(tsv_file)
    strain_name = filename.replace("_blast.tsv", "")
    out_var = os.path.join(out_dir, f"{strain_name}.var")
    out_stops = os.path.join(out_dir, f"{strain_name}.stops")
    stats = Counter()

    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    try:
        # 1) Lecture des hits
        hits = []
        with open(tsv_file, "r") as f_tsv:
            for line in f_tsv:
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) < 8:
                    stats["skipped_malformed"] += 1
                    continue
                hits.append(parts)

        # 2) Statistiques de collisions (plusieurs CDS -> même gène de référence)
        stats["n_hits"] = len(hits)
        stats["n_ref_genes"] = len({p[1] for p in hits})
        stats["collisions"] = stats["n_hits"] - stats["n_ref_genes"]

        # 3) Option : un seul CDS par gène de référence (meilleur bitscore)
        if best_per_gene:
            best = {}
            for p in hits:
                if p[1] not in best or float(p[5]) > float(best[p[1]][5]):
                    best[p[1]] = p
            stats["dropped_collisions"] = len(hits) - len(best)
            hits = list(best.values())

        # 4) Extraction
        var_lines, stop_lines = [], []
        for p in hits:
            qid, sseqid, qseq, sseq = p[0], p[1], p[6], p[7]
            try:
                sstart = int(p[8]) if len(p) > 8 else 1
            except ValueError:
                stats["skipped_malformed"] += 1
                continue

            if qseq == sseq:
                continue                                   # identique : aucun variant

            if len(sseq) != len(qseq) and '-' not in sseq and '-' not in qseq:
                # champs DIAMOND sans gaps : on retrouve les indels par réalignement
                aln = aligner.align(sseq, qseq)[0]
                sseq, qseq = str(aln[0]), str(aln[1])
                stats["realigned"] += 1

            if len(sseq) != len(qseq):                     # cas incohérent : on ne devine pas
                stats["skipped_length_mismatch"] += 1
                continue

            muts, stops, st = extract_mutations(sseq, qseq, ref_start=sstart,
                                                truncate_at_stop=truncate_at_stop)
            stats.update(st)

            suffix = f"\t{qid}" if with_qseqid else ""
            var_lines.extend(f"{sseqid}\t{m}{suffix}" for m in muts)
            stop_lines.extend(f"{sseqid}\t{s}{suffix}" for s in stops)

        stats["n_mutations"] = len(var_lines)

        with open(out_var, "w") as f_out:
            if var_lines:
                f_out.write("\n".join(var_lines) + "\n")
        with open(out_stops, "w") as f_stops:
            if stop_lines:
                f_stops.write("\n".join(stop_lines) + "\n")

        return strain_name, dict(stats), None
    except Exception as e:
        return strain_name, dict(stats), str(e)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extracts mutations and indels from TSV alignments.")

    parser.add_argument("-t", "--tsv_dir", default="/data/projet1/02_Pipeline_Output/2_Blast")
    parser.add_argument("-o", "--out_dir", default="/data/projet1/02_Pipeline_Output/3_Mutations")
    parser.add_argument("-j", "--jobs", type=int, default=4)
    parser.add_argument("--best-per-gene", action="store_true",
                        help="Ne garder qu'un hit (meilleur bitscore) par gène de référence et par fichier.")
    parser.add_argument("--with-qseqid", action="store_true",
                        help="Ajoute le qseqid en 3e colonne des .var et .stops.")
    parser.add_argument("--truncate-at-stop", action="store_true",
                        help="Arrête l'extraction d'un CDS au premier stop prématuré.")

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    tsv_files = sorted(glob.glob(os.path.join(args.tsv_dir, "*.tsv")))

    if not tsv_files:
        print(f"Error: No TSV files found in {args.tsv_dir}", flush=True)
        exit(1)

    print(f"Running extraction on {len(tsv_files)} files with {args.jobs} cores...", flush=True)

    tasks = [(tsv, args.out_dir, args.best_per_gene, args.with_qseqid, args.truncate_at_stop)
             for tsv in tsv_files]

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

    print(f"\nStep 3 completed: {completed} .var files generated in {args.out_dir} "
          f"({total['n_mutations']} mutations total; {total['stop_gain']} premature stops "
          f"written to .stops files; {errors} errors).", flush=True)
    print("Detail: " + ", ".join(f"{c}={total[c]}" for c in STAT_COLS if c != "n_mutations"), flush=True)
    print(f"Per-file statistics: {stats_path}", flush=True)