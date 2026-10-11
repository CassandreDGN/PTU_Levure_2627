# Projet PTU Bioinformatique — Analyse des Variantes et Prédiction Fonctionnelle

## 1. Description et Objectif du Projet
L'objectif du projet est d'analyser les CDS et séquences protéiques chez *Saccharomyces cerevisiae*, d'effectuer des alignements de séquences (BLASTn ou p) et d'évaluer l'impact fonctionnel des mutations identifiées à l'aide de l'outil PROVEAN.

---

## 2. Organisation du workspace

Voici les explications de l'organisation du workspace : 

```text
.
├── 00_Data/                # Données brutes et bases de données (CDS, références, bases PROVEAN)
├── 01_Script/              # Scripts Python et bash de traitement (traduction, BLAST, extraction de mutations)
├── 02_Pipeline_Output/     # Résultats intermédiaires des différentes étapes  et logs d'exécution (tous les outputs des scripts en gros)
├── 03_Archive/             # Fichiers archivés, anciennes versions, de scripts et tests, fichiers ne trouvant par leur place ailelurs
├── 04_Resultats_BlastN/    # Résultats spécifiques aux alignements BLASTn (branche expérimentale de notre projet mené par Beray )
├── 05_Documentation/       # Cahier de laboratoire et fichiers d'environnement Conda (.yml) par exemple
├── 06_Conda/               # Environnements virtuels locaux (provean_env, env_blast)
└── README.md               # Présentation générale du projet et guides d'exécution