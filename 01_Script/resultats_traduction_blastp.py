from Bio import SeqIO
import os

# Dossier source et dossier de sortie isolé
dossier_entree = "00_Data/cds_fasta/"
dossier_sortie = "resultats_proteines/"

os.makedirs(dossier_sortie, exist_ok=True)

print("--- Etape 1 : Traduction CDS -> Proteines & Controle qualite ---")

total_seq = 0
filtrees_taille = 0
stops_internes = 0

# Parcours des fichiers fasta / fna dans le dossier source
for filename in os.listdir(dossier_entree):
    if filename.endswith((".fasta", ".fa", ".fna")):
        chemin_in = os.path.join(dossier_entree, filename)
        chemin_out = os.path.join(dossier_sortie, "prot_" + filename)
        
        with open(chemin_out, "w") as out_f:
            for record in SeqIO.parse(chemin_in, "fasta"):
                total_seq += 1
                
                # Traduction complète (to_stop=False pour garder les stops et les repérer)
                prot_seq = record.seq.translate(to_stop=False)
                
                # Contrôle : détection d'un codon stop interne (un '*' avant la fin)
                if '*' in prot_seq[:-1]:
                    stops_internes += 1
                
                # Contrôle taille : élimination des protéines anormalement courtes (< 30 AA)
                if len(prot_seq) >= 30:
                    record.seq = prot_seq
                    SeqIO.write(record, out_f, "fasta")
                else:
                    filtrees_taille += 1

print("--- Etape 1 terminee ---")
print(f"Total sequences analysees : {total_seq}")
print(f"Rejetees (trop courtes < 30 AA) : {filtrees_taille}")
print(f"Proteines avec stop(s) interne(s) : {stops_internes}")
print(f"Resultats propres stockes dans : {dossier_sortie}")
