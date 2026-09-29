# Cahier de laboratoire - PTU Levures

## 13/09/2026

- Mise en place du GitHub par Cassandre

***

## 15/09/2026

* **Beray** : Préparation de l'espace de travail, téléchargement et indexation du génome de référence S288C, réalisation du test unitaire et exécution du BLAST sur l'ensemble des souches.
* **Cassandre** : Documentation pour reproductibilité

>### Environnement logiciel & Outils
* **Système** : Serveur Linux (`/data/projet1`)
* **Gestionnaire d'environnement** : Conda
* **Outils principaux** :
  * `blastn` (NCBI BLAST+ v2.x via environnement Conda `env_blast`)
  * `makeblastdb`


>### Lignes de commande utilisées

### Etape 1 : Mise en place
* Répertoire de travail principal : `/data/projet1/ptu_projet`
* Données d'entrée : 400 souches diploïdes de *Saccharomyces cerevisiae* phasées (2 fichiers `.fa` par souche).
* Les identifiants de gènes bruts nécessitant une étape d'annotation par alignement.

### Etape 2 : Téléchargement du génome de référence
Création d'un répertoire mutualisé et récupération des CDS officielles de *S. cerevisiae* depuis SGD :

```bash
cd /data/projet1
mkdir -p reference_S288C
cd reference_S288C

# Téléchargement et décompression de la référence
wget http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding_all.fasta.gz
gunzip orf_coding_all.fasta.gz
```

### Etape 3 : Configuration de l'environnement logiciel et indexation de la base

* Création et activation de l'environnement
```bash
conda create -n env_blast -c bioconda blast -y
conda activate env_blast
```
* Indexation de la base de données nucléotidique
```bash
makeblastdb -in orf_coding_all.fasta -dbtype nucl -title "Ref_S288C" -out Ref_S288C
```
### Etape 4 : Test unitaire d'annotation BLASTN (validation du paramétrage sur un fichier test)

* Identification d'un fichier cobaye
```bash
ls *.fa | head -n 5
```
* Lancement du test unitaire
```bash
blastn -query AAC.HP1.

nuclear_genome.Final.cds.fa \
       -db /data/projet1/reference_S288C/Ref_S288C \
       -out resultat_test.txt \
       -outfmt 6 \
       -max_target_seqs 1
```

* Vérification du fichier de sortie (~621 Ko) et inspection
```bash
ls -lh resultat_test.txt
head -n 5 resultat_test.txt
```

### Phase 5 : Automatisation de l'annotation à pour toutes les données

```bash
mkdir -p /data/projet1/resultats_blast
cd /data/projet1/ptu_projet

for fichier in *.fa; do
    nom_base=$(basename "$fichier" .fa)
    blastn -query "$fichier" \
           -db /data/projet1/reference_S288C/Ref_S288C \
           -out "/data/projet1/resultats_blast/${nom_base}_resultat.txt" \
           -outfmt 6 \
           -max_target_seqs 1
    echo "Analyse terminée pour : $fichier"
done
```

***

## 28/09/2026
* **Rayan** : Ecriture d'un plan déterminant les étapes à suivre pour le déroulement du projet
* **Cassandre** : Développement du script de traduction des CDS en protéines (Étape 1 du pipeline) et mise en place du contrôle qualité (QC) strict sur un sous-échantillon de test.

>### Environnement logiciel & Outils
* **Langage** : Python
* **Outils principaux** : 
  * Script `1_CDStoPROTEIN.py` (librairies standards `os`, `sys`, `argparse`)

>### Lignes de commande utilisées

### Etape 1 : Développement du script de traduction et de QC (`1_CDStoPROTEIN.py`)
* Création d'un script Python permettant la traduction nucléotide vers acide aminé.
* Implémentation des règles de rejet strictes (QC) pour isoler les protéines valides :
  * Longueur non multiple de 3 (frameshift)
  * Absence de codon START (`ATG`)
  * Présence de nucléotides ambigus (`N`) ou codons invalides
  * Présence d'un codon STOP prématuré
  * Protéines trop courtes (< 30 AA)

### Etape 2 : Test sur un sous-échantillon
* Données d'entrée (test) : `/data/projet1/subgroup_test` (6 fichiers `.cds.fa`).
* Lancement du script sur les données de test :

```bash
python3 /data/projet1/01_Script/1_CDStoPROTEIN.py \
  -i /data/projet1/subgroup_test \
  -o /data/projet1/02_Test_Proteins \
  -m 30
```

### Etape 3 : Débogage et Troubleshooting
* **Problème rencontré** : Taux de rejet >85% avec le motif `Codon inv`. 
* **Inspection des séquences rejetées** :

```bash
head -n 11 /data/projet1/02_Test_Proteins/rejected_sequences.tsv
```

* **Diagnostic et Correction** : Le log d'erreur a montré que le codon `TGG` était systématiquement rejeté. Le dictionnaire interne `CODON_TABLE` a été mis à jour pour intégrer le Tryptophane (W) manquant et sécurisé avec les 64 codons complets.

### Etape 4 : Exécution finale et validation
* Relance du script après modification du dictionnaire.
* Vérification des sorties :
  * Génération des séquences protéiques propres (`*.prot.fasta`).
  * Création des journaux de logs : `qc_summary.tsv` (confirmant le passage de la quasi-totalité des CDS) et `rejected_sequences.tsv` (détaillant la position des quelques anomalies restantes).
