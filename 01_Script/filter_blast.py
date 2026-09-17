# filter_blast.py
import argparse
import os


def parse__blast(input_file, output_file):
    best_matches = {}

    with open(input_file, "r") as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) < 2:
                continue

            gene_souche = cols[0].split("|")[0]
            gene_ref = cols[1]

            best_matches[gene_souche] = gene_ref #peut être juste ignorer si clé déjà existante

    with open(output_file, "w") as out:
        out.write("gene_souche\tgene_ref\n")
        for g_souche, g_ref in best_matches.items():
            out.write(f"{g_souche}\t{g_ref}\n")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Parcourt un dossier de résultats BLAST et génère un fichier TSV filtré pour chaque souche."
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Dossier contenant les fichiers de résultats BLAST en input",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Dossier de destination des fichiers TSV",
    )

    args = parser.parse_args()

    input_dir = args.input
    output_dir = args.output if args.output else input_dir

    os.makedirs(output_dir, exist_ok=True)

    for filename in os.listdir(input_dir):
        if filename.endswith(".txt"):
            input_path = os.path.join(input_dir, filename)

            strain_prefix = ".".join(filename.split(".")[:2])
            generated_filename = f"clean_results_{strain_prefix}.tsv"
            output_path = os.path.join(output_dir, generated_filename)

            print(f"Traitement de {filename} -> {generated_filename}...")
            parse__blast(input_path, output_path)

    print("Fin du traitement.")