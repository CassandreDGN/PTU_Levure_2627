# Cahier de laboratoire - PTU Levures
# Versioning : https://github.com/CassandreDGN/PTU_Levure_2627

## 13/09/2026

- Mise en place du GitHub (https://github.com/CassandreDGN/PTU_Levure_2627) par Cassandre

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

- **Rayan** : Écriture du plan déterminant les étapes à suivre pour le projet.
- **Cassandre** : Développement du script de traduction des séquences (CDS vers protéines) et mise en place du contrôle qualité (QC).

>### Environnement logiciel & Outils
- **Langage** : Python 3.14.7
- **Script créé** : `01_Script/1_CDStoPROTEIN.py` (librairies `os`, `sys`, `argparse`)

>### Travail effectué et Lignes de commande

- **Développement des critères QC** : Implémentation des règles de rejet strictes (frameshift, absence de codon START, nucléotides/codons invalides, STOP prématuré, ou séquence < 30 AA).
- **Débogage sur sous-échantillon** : Test sur 6 fichiers cibles. Identification et correction d'un taux de rejet anormal (>85%) lié à l'omission du codon Tryptophane (`TGG`) dans le dictionnaire interne du script.
- **Validation** : Relance du script corrigé, génération des séquences protéiques propres (`*.prot.fasta`) et des fichiers de suivi (`qc_summary.tsv`, `rejected_sequences.tsv`).

```bash
python3 /data/projet1/01_Script/1_CDStoPROTEIN.py \
  -i /data/projet1/subgroup_test \
  -o /data/projet1/02_Pipline_Output \
  -m 30
```

## 06/10/2026

- **Rayan & Cassandre** : Développement du script d'alignement BLAST pour les séquences protéiques.
- **Cassandre** : Création de l'environnement logiciel dédié.

>### Environnement logiciel & Outils
- **Gestionnaire** : Conda (`my_env_blast` incluant Biopython)
- **Script créé** : `01_Script/02.1_blastP.py` (librairies `os`, `glob`, `argparse`, `subprocess`)

>### Travail effectué et Lignes de commande

- **Automatisation BLASTp** : Ébauche par Rayan et finalisation par Cassandre d'un script exécutant `blastp` pour chaque souche contre la base de référence S288C.
- **Génération des résultats** :  formatage tabulaire personnalisé pour extraire les alignements significatifs.

```bash
python3 /data/projet1/01_Script/02.1_blastP.py
```
## 07/10/2026

- **Cassandre** : Développement du script d'extraction des mutations au format PROVEAN (`.var`) à partir des alignements BLASTp.
- **Beray** : Traduction des séquences pour adaptation au cheminement et validation du flux BLASTn.

### Script d'extraction des mutations
>### Environnement logiciel & Outils
- **Gestionnaire** : Conda (`blastp_env`)
- **Langage et Librairies** : Python 3.14.7 (Biopython : `SeqIO`, `Align.PairwiseAligner`)
- **Script créé** : `01_Script/03.1_extract_mutations.py`

>### Travail effectué et Lignes de commande

- **Extraction des variants** : Utilisation de `PairwiseAligner` (mode global) pour aligner les séquences protéiques des souches contre la référence S288C (`orf_trans_all.fasta`).
- **Formatage PROVEAN** : Detection et classification automatique des variations (substitutions, insertions, délétions) générant un fichier `.var` par souche.
- **Exécution du pipeline** :
  ```bash
  python3 /data/projet1/01_Script/03.1_extract_mutations.py \
    -t /data/projet1/02_Pipeline_Output/02_Blast \
    -q /data/projet1/02_Pipeline_Output/traduction_test_batch \
    -r /data/projet1/00_Data/reference_S288C/orf_trans_all.fasta \
    -o /data/projet1/02_Pipeline_Output/03_Mutations

##10 octobre 2026
-  **Cassandre** 

### Objectif
Exécution et automatisation du pipeline d'analyse fonctionnelle des variantes protéiques avec PROVEAN.

### Travail effectué & Problèmes rencontrés
- **Tentatives d'exécution du script PROVEAN :**
  - Configuration et tests prolongés du script de lancement PROVEAN (`03.2_run_provean.sh`) et des dépendances associées (PSI-BLAST, bases de données de référence).
  - Nombreux échecs et erreurs rencontrés lors des exécutions successives sur le jeu de données, principalement liés à la gestion des chemins, des formats d'entrée et des dépendances logicielles.
  - Après plusieurs heures de sueur et de rage intense sur les multiples essais, le ragequit a gagné et la  décision a été prise de nettoyer le répertoire de travail et d'éliminer le script dysfonctionnel (`03.2_run_provean.sh`) afin de maintenir un espace de travail propre et exempt de scripts obsolètes.

### État de l'environnement & Perspectives
- **Conservation de l'environnement virtuel :**
  - L'environnement Conda dédié (`provean_env`), situé sous `06_Conda/provean_env`, a été conservé intact sur la machine virtuelle.