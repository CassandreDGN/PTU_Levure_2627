#!/usr/bin/env python3

import os
import glob
import argparse
from Bio import Align
from concurrent.futures import ProcessPoolExecutor, as_completed

def extract_mutations(sseq, qseq, ref_start=1):
    mutations = []
    ref_pos = ref_start - 1
    
    for r, q in zip(sseq, qseq):
        if r != '-':
            ref_pos += 1
        
        if r in ('*', 'X') or q in ('*', 'X'):
            continue
            
        if r != q:
            if r == '-':  # Insertion
                mutations.append(f"{ref_pos}_{ref_pos+1}ins{q}")
            elif q == '-':  # Deletion
                mutations.append(f"{r}{ref_pos}del")
            else:  # Substitution
                mutations.append(f"{r}{ref_pos}{q}")
                
    return mutations

def process_tsv_single(args_tuple):
    tsv_file, out_dir = args_tuple
    filename = os.path.basename(tsv_file)
    strain_name = filename.replace("_blast.tsv", "")
    out_var = os.path.join(out_dir, f"{strain_name}.var")

    # Local aligner instance for fallback handling of indels
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    mutations_count = 0

    try:
        with open(tsv_file, 'r') as f_tsv, open(out_var, 'w') as f_out:
            for line in f_tsv:
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                
                sseqid = parts[1]
                qseq = parts[6]
                sseq = parts[7]
                sstart = int(parts[8]) if len(parts) > 8 else 1

                if sseq == qseq:
                    continue  # Identical, no mutations

                # If lengths differ and no gaps exist, realign snippet with Biopython to insert '-'
                if len(sseq) != len(qseq) and '-' not in sseq and '-' not in qseq:
                    alignments = aligner.align(sseq, qseq)
                    best_aln = alignments[0]
                    sseq_aln = str(best_aln[0])
                    qseq_aln = str(best_aln[1])
                else:
                    sseq_aln = sseq
                    qseq_aln = qseq

                mutations = extract_mutations(sseq_aln, qseq_aln, ref_start=sstart)

                for mutant in mutations:
                    f_out.write(f"{sseqid}\t{mutant}\n")
                    mutations_count += 1

        return strain_name, mutations_count, None
    except Exception as e:
        return strain_name, 0, str(e)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extracts mutations and indels from TSV alignments.")

    parser.add_argument("-t", "--tsv_dir", default="/data/projet1/02_Pipeline_Output/2_Blast")
    parser.add_argument("-o", "--out_dir", default="/data/projet1/02_Pipeline_Output/3_Mutations")
    parser.add_argument("-j", "--jobs", type=int, default=4)

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    tsv_files = sorted(glob.glob(os.path.join(args.tsv_dir, "*.tsv")))

    if not tsv_files:
        print(f"Error: No TSV files found in {args.tsv_dir}", flush=True)
        exit(1)

    print(f"Running extraction on {len(tsv_files)} files with {args.jobs} cores...", flush=True)

    tasks = [(tsv, args.out_dir) for tsv in tsv_files]

    completed = 0
    total_mutations = 0
    errors = 0

    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(process_tsv_single, task) for task in tasks]
        for future in as_completed(futures):
            strain_name, count, err = future.result()
            if err:
                errors += 1
                print(f"[{completed + errors}/{len(tsv_files)}] {strain_name} -> ERROR: {err}", flush=True)
            else:
                completed += 1
                total_mutations += count

            if (completed + errors) % 50 == 0 or (completed + errors) == len(tsv_files):
                print(f"Progress: {completed + errors}/{len(tsv_files)} files processed...", flush=True)

    print(f"\nStep 3 completed: {completed} .var files generated in {args.out_dir} ({total_mutations} mutations total).", flush=True)