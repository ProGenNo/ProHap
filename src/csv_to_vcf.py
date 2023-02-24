import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Formats a CSV file to fit the VCF requirements")

parser.add_argument("-i", dest="input_file", required=True,
                    help="input file", metavar="FILE")

parser.add_argument("-sep", dest="sep", required=False, type=str,
                    help="separator in the input file", default=',')

parser.add_argument("-ch", dest="chrom_column", required=False,
                    help="chromosome column header", default="chrom")

parser.add_argument("-pos", dest="pos_column", required=False,
                    help="position column header", default="pos")

parser.add_argument("-id", dest="id_column", required=False,
                    help="ID column header", default="accession")

parser.add_argument("-r", dest="ref_column", required=False,
                    help="ref. allele column header", default="ref")

parser.add_argument("-a", dest="alt_column", required=False,
                    help="alt.allele column header", default="alt")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output VCF file", metavar="FILE")

args = parser.parse_args()

df = pd.read_csv(args.input_file, sep=args.sep)

df['#CHROM'] = df[args.chrom_column]
df['POS'] = df[args.pos_column]
df['ID'] = df[args.id_column]
df['REF'] = df[args.ref_column]
df['ALT'] = df[args.alt_column]
df['INFO'] = '.'
df['QUAL'] = '.'
df['FILTER'] = '.'

df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']].sort_values(by=['#CHROM', 'POS']).to_csv(args.output_file, sep='\t', header=True, index=False)