import pandas as pd
import argparse
import re

parser = argparse.ArgumentParser(
    description='Makes a list of variants from the annotated list of peptides')

parser.add_argument("-i", dest="input_filename", required=True,
                    help="input file")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output CSV file", metavar="FILE")

args = parser.parse_args()

df = pd.read_csv(args.input_filename)

variants_data = {}
for index,row in df.iterrows():
    if (row['covered_changes_dna'] == row['covered_changes_dna']):
            for change in re.split(r"[;|]", row['covered_changes_dna']):
                    if change in variants_data:
                            variants_data[change] += 1
                    else:
                            variants_data[change] = 1

df_data = [ [change, variants_data[change]] for change in list(variants_data.keys()) ]
variants_df = pd.DataFrame(data=df_data, columns=['Variant', 'matching_peptides'])
variants_df.to_csv(args.output_file, header=True, index=False)
