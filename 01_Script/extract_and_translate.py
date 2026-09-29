"""convertir les cds en séquences protéiques
    relier les 2 dans un dictionnaire de dictionnaire
    ex 

    protein_data = {
    "YHL050C": {
        "gene_souche": "AAC_00012",
        "seq_souche": "MSLVTK...",
        "seq_ref": "MSLVTK...",
    },
    "YGR192C": {
        "gene_souche": "AAC_00045",
        "seq_souche": "MAAAGT...",
        "seq_ref": "MAAGTT...",
    },
}

à partir de ça, write dans un fichier tsv pour utiliser en input de provean
"""
import argparse
import os 