from Bio import SeqIO
import os

dossier_entree = "Data/"
dossier_sortie = "02_Test_Proteins/"

os.makedirs(dossier_sortie, exist_ok=True)

print("Début de la traduction...")

for filename in os.listdir(dossier_entree):
  if filename.endswith(".fasta") or filename.e          SeqIO.write(record, out_f, "fasta")

print("Traduction terminée !")
