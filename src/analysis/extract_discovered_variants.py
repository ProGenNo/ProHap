import pandas as pd
import argparse
import re

parser = argparse.ArgumentParser(
    description='Makes a list of variants from the annotated list of peptides')

parser.add_argument("-i", dest="input_filename", required=True,
                    help="input file")

parser.add_argument("-hap", dest="haplo_db", required=False,
                    help="haplotypes tab-separated file (optional)", default=None)
                    
parser.add_argument("-hap_prefix", dest="haplo_prefix", required=False,
                    help="prefix for haplotype protein ID (default: 'haplo_')", default='haplo_')

parser.add_argument("-o", dest="output_file", required=True,
                    help="output CSV file", metavar="FILE")

args = parser.parse_args()

df = pd.read_csv(args.input_filename)

haplo_db = None
if (args.haplo_db):
    print ("Reading", args.haplo_db)
    haplo_db = pd.read_csv(args.haplo_db, sep='\t', header=0)
    haplo_db.set_index('HaplotypeID', inplace=True)

variants_data = {}
for index,row in df.iterrows():
    if (row['covered_changes_dna'] == row['covered_changes_dna']):

            max_freq = -1

            if (haplo_db is not None):
                max_freq = max([ haplo_db.loc[haploID]['frequency'] for haploID in row['matching_proteins'].split(';') if haploID.startswith(args.haplo_prefix) ])
            
            for change in re.split(r"[;|]", row['covered_changes_dna']):
                    if change in variants_data:
                            variants_data[change]['count'] += 1
                            variants_data[change]['max_freq'] = max(variants_data[change]['max_freq'], max_freq)
                    else:
                            variants_data[change] = { 'count': 1, 'max_freq': max_freq }

df_data = [ [change, variants_data[change]['count'], variants_data[change]['max_freq']] for change in list(variants_data.keys()) ]
variants_df = pd.DataFrame(data=df_data, columns=['Variant', 'matching_peptides', 'max_haplotype_frequency'])
variants_df.to_csv(args.output_file, header=True, index=False)
