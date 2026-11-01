#!/usr/bin/env python3
import os
import concurrent.futures
from pathlib import Path

MUTATIONS_DIR = Path("/data/projet1/02_Pipeline_Output/3_SIFT_Mutations")
SIFT_RESULTS_DIR = Path("/data/projet1/02_Pipeline_Output/4_SIFT_SWISSPROT_Results")

def process_haplotype(hap_dir):
    hap_name = hap_dir.name
    base_strain = hap_name.split(".nuclear_genome")[0]
    
    sift_hap_dir = SIFT_RESULTS_DIR / hap_name
    if not sift_hap_dir.exists():
        sift_hap_dir = SIFT_RESULTS_DIR / base_strain

    if not sift_hap_dir.exists():
        return f"[SKIP] {hap_name} (Aucun dossier)"

    subst_hap_dir = hap_dir / "subst" 

    # 1. Événements LoF
    lof_info = {}
    lof_file = hap_dir / "loss_of_function_events.tsv"
    if lof_file.exists():
        with open(lof_file, "r") as f_lof:
            lines = f_lof.read().splitlines()
            if len(lines) > 1:
                for line in lines[1:]:
                    parts = line.split("\t")
                    if len(parts) >= 2:
                        lof_info[parts[0]] = parts[1]

    # 2. Prédictions SIFT4G
    gene_summary = {}
    
    # OPTIMISATION 1 : os.scandir est immensément plus rapide que pathlib.glob
    try:
        sift_files = [f for f in os.scandir(sift_hap_dir) if f.name.endswith(".SIFTprediction")]
    except FileNotFoundError:
        return f"[SKIP] {hap_name}"

    for entry in sift_files:
        prot_id = entry.name.replace(".SIFTprediction", "")
        
        # OPTIMISATION 2 : Lecture du bloc complet et méthode count() en C (ultra-rapide)
        with open(entry.path, "r") as f_sift:
            content = f_sift.read()
            n_del = content.count("DELETERIOUS")
            n_tol = content.count("TOLERATED")
            n_subst = n_del + n_tol
        
        has_stop = "YES" if prot_id in lof_info and "STOP" in lof_info[prot_id] else "NO"
        has_large_del = "YES" if prot_id in lof_info and "DELETION" in lof_info[prot_id] else "NO"
        
        if prot_id in lof_info:
            is_functional = "NO"
            detailed_status = "LOF_STOP" if has_stop == "YES" else "LOF_DELETION"
        elif n_del > 0:
            is_functional = "NO"
            detailed_status = "DELETERIOUS_SUBST"
        elif n_tol > 0 and n_del == 0:
            is_functional = "YES"
            detailed_status = "TOLERATED_SUBST"
        elif n_subst == 0:
            subst_file = subst_hap_dir / f"{prot_id}.subst"
            if subst_file.exists() and subst_file.stat().st_size > 0:
                is_functional = "UNKNOWN"
                detailed_status = "SIFT_FAILED"
            else:
                is_functional = "YES"
                detailed_status = "INTACT_WILDTYPE"
        else:
            is_functional = "UNKNOWN"
            detailed_status = "UNKNOWN"
            
        gene_summary[prot_id] = {
            "n_subst": n_subst, "n_tol": n_tol, "n_del": n_del,
            "has_stop": has_stop, "has_large_del": has_large_del,
            "detailed_status": detailed_status, "is_functional": is_functional
        }

    for prot_id, lof_type in lof_info.items():
        if prot_id not in gene_summary:
            has_stop = "YES" if "STOP" in lof_type else "NO"
            has_large_del = "YES" if "DELETION" in lof_type else "NO"
            gene_summary[prot_id] = {
                "n_subst": 0, "n_tol": 0, "n_del": 0,
                "has_stop": has_stop, "has_large_del": has_large_del,
                "detailed_status": "LOF_STOP" if has_stop == "YES" else "LOF_DELETION",
                "is_functional": "NO"
            }

    # 3. Écriture en une seule passe
    summary_file = hap_dir / "gene_summary_swissprot.tsv"
    out_lines = ["protein\tn_subst\tn_tolerated\tn_deleterious\thas_stop\thas_large_del\tdetailed_status\tis_functional\n"]
    for prot_id, data in sorted(gene_summary.items()):
        out_lines.append(f"{prot_id}\t{data['n_subst']}\t{data['n_tol']}\t{data['n_del']}\t{data['has_stop']}\t{data['has_large_del']}\t{data['detailed_status']}\t{data['is_functional']}\n")
    
    with open(summary_file, "w") as fout:
        fout.writelines(out_lines)

    return f"[OK] {hap_name} ({len(gene_summary)} gènes)"

if __name__ == "__main__":
    print("=== ÉTAPE 3.5 ULTRA-RAPIDE : AGRÉGATION SIFT4G & LOF ===")
    haplotypes = [d for d in MUTATIONS_DIR.iterdir() if d.is_dir()]
    
    # OPTIMISATION 3 : Affichage as_completed pour bypasser le buffer
    with concurrent.futures.ProcessPoolExecutor(max_workers=3) as executor:
        futures = {executor.submit(process_haplotype, hap): hap for hap in haplotypes}
        for future in concurrent.futures.as_completed(futures):
            res = future.result()
            if "[OK]" in res:
                print(res, flush=True)