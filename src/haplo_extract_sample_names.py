import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Reads a TSV file produced by ProHap, extracts the "samples" column to a separate file and removes it from the original')

parser.add_argument("-hap_tsv", dest="haplo_table", required=True,
                    help="input haplotype table", metavar="FILE")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output TSV file", metavar="FILE")

parser.add_argument("-samples", dest="output_file_samples", required=True,
                    help="output file (tab-separated) for the list of samples per haplotype", metavar="FILE")

args = parser.parse_args()

haplo_df = pd.read_table(args.haplo_table, compression='infer')

haplo_df.drop('samples', axis=1).to_csv(args.output_file, sep='\t', index=False, compression='infer')
haplo_df[['HaplotypeID', 'samples']].to_csv(args.output_file_samples, sep='\t', index=False, compression='infer')