#!/usr/bin/env python3
"""
Classification sens / non-sens des variants (fichiers .var du projet PTU).
Lecture seule : aucun fichier n'est modifie, ni supprime, ni deplace.

Format des fichiers .var (deduit des donnees observees sur le serveur ptu.bigest-icube.fr) :
  - substitution simple       : AA + position + AA  (ex: Q487K, L63S, Y351N)
  - deletion                  : AA + position + del (ex: F177del, T655del)
  - insertion                 : pos1_pos2insAA    (ex: 77_78insD)
  - codon stop (nonsense)     : *                   (ex: Q487*, W23*)
"""
import sys, time, re
import os

HOST = "ptu.bigest-icube.fr"
USER = "sonmez"
PWD  = "yfmPFKczpX"

import paramiko
client = paramiko.SSHClient()
client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
client.connect(HOST, username=USER, password=PWD, timeout=20,
               allow_agent=False, look_for_keys=False)

import os
def run(cmd, timeout=1800):
    stdin, stdout, stderr = client.exec_command(cmd, timeout=timeout)
    out = stdout.read().decode(errors="replace")
    err = stderr.read().decode(errors="replace")
    rc = stdout.channel.recv_exit_status()
    return out, err, rc

ROOT = "/data/projet1/02_Pipeline_Output/3_Mutations"

cmd_list = f"find {ROOT} -type f -name '*.var' | sort"
out, err, rc = run(cmd_list)
if rc != 0:
    print("ERREUR liste fichiers:", err[:300], file=sys.stderr)
    sys.exit(1)
files = [l for l in out.splitlines() if l.strip()]
print(f"Nombre de fichiers .var trouves: {len(files)}", flush=True)

SUB_RE    = re.compile(r"^[A-Z][0-9]+[A-Z]$")
DEL_RE    = re.compile(r"^[A-Z][0-9]+del$")
INS_RE    = re.compile(r"^[0-9]+_[0-9]+ins[A-Z]+$")

counts = {"substitution": 0, "deletion": 0, "insertion": 0, "nonsense": 0, "other": 0}
nonsens_examples = []
ex_sub = []; ex_del=[]; ex_ins=[]

start = time.time()
for idx, f in enumerate(files):
    souche = f.split("/")[-1].replace(".var", "")
    cmd_cat = f"cat '{f}'"
    out, err, rc = run(cmd_cat)
    if rc != 0:
        print(f"  ERREUR lecture {souche}: {err[:200]}", flush=True)
        continue
    for line in out.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        var = parts[1].strip()
        if SUB_RE.match(var):
            counts["substitution"] += 1
            ex_sub.append((souche, var))
        elif DEL_RE.match(var):
            counts["deletion"] += 1
            ex_del.append((souche, var))
        elif INS_RE.match(var):
            counts["insertion"] += 1
            ex_ins.append((souche, var))
        elif "*" in var:
            counts["nonsense"] += 1
            nonsens_examples.append((souche, var))
        else:
            counts["other"] += 1
    if (idx + 1) % 50 == 0:
        print(f"  {idx+1}/{len(files)} fichiers, elapsed {time.time()-start:.0f}s",
              flush=True)

total = sum(counts.values())
print("\n=== RESULTATS FINAUX ===", flush=True)
print("totaux:", counts, flush=True)
print("TOTAL:", total, flush=True)
if total > 0:
    for k, v in counts.items():
        print(f"  {k}: {v} ({100*v/total:.2f}%)", flush=True)
if nonsens_examples:
    print("exemples nonsense:", nonsens_examples[:10], flush=True)
print("ex_sub sample:", ex_sub[:5], flush=True)
print("ex_del sample:", ex_del[:5], flush=True)
print("ex_ins sample:", ex_ins[:5], flush=True)

client.close()
