# filter_blast.py
import argparse
import os


def parse__blast(input_file, output_file):
    best_matches = {}

    with open(input_file, "r") as f:
        for line in f:
            line = line.strip()


            cols = line.split("\t")
            if len(cols) < 12:
                continue

            gene_souche = cols[0].split("|")[0]
            gene_ref = cols[1]
            bit_score = float(cols[11])

            if (gene_souche not in best_matches or bit_score > best_matches[gene_souche][1]):best_matches[gene_souche] = (gene_ref, bit_score)

    with open(output_file, "w") as out:
        out.write("gene_souche\tgene_ref\n")
        for g_souche, (g_ref, x) in best_matches.items():
            out.write(f"{g_souche}\t{g_ref}\n")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Nettoie un fichier output de Blast pour retenir seulement le meilleur match souche/reference pour chaque gène ainsi que le score")

    parser.add_argument("-i", "--input", required=True, help="Fichier .txt résultat BLAST en input")
    parser.add_argument("-o", "--output", help="Fichier TSV de sortie ou dossier de destination")

    args = parser.parse_args()

    base_filename = os.path.basename(args.input)
    strain_prefix = ".".join(base_filename.split(".")[:2])
    generated_filename = f"mapping_{strain_prefix}.tsv"

    if not args.output:
        args.output = generated_filename
    elif os.path.isdir(args.output):
        args.output = os.path.join(args.output, generated_filename)

    parse__blast(args.input, args.output)