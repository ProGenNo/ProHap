import argparse
import pandas as pd
import re

parser = argparse.ArgumentParser(description="Splits a VCF into separate files for each chromosome, fills in missing values by '-'.")

parser.add_argument("-i", dest="input_file", required=True,
                    help="input CSV file", metavar="FILE")

parser.add_argument("-h", dest="header_size", required=False, type=int,
                    help="number of lines of the VCF header", default=1)

parser.add_argument("-sep", dest="infile_sep", required=False, type=str,
                    help="separator in the VCF file", default=',')

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

parser.add_argument("-o", dest="output_file_prefix", required=True,
                    help="output VCF file prefix", metavar="FILE")

args = parser.parse_args()

# if (not re.match(r'[CGTA]', str(ALT)))

df = pd.read_csv(args.input_file, sep=args.infile_sep, header=args.header_size-1)
CHROMOSOMES = [str(x) for x in list(range(1, 23))] + ['X', 'Y']

def fill_missing(val):
	if (not re.match(r'[CGTA]', str(val))):
		return '-'
	return val

df['#CHROM'] = df[args.chrom_column].apply(lambda x: str(x).replace('chr', ''))
df['POS'] = df[args.pos_column]
df['ID'] = df[args.id_column]
df['REF'] = df[args.ref_column].apply(fill_missing)
df['ALT'] = df[args.alt_column].apply(fill_missing)

df[['#CHROM','POS','ID','REF','ALT']].to_csv(args.output_file_prefix + '.vcf', sep='\t', index=False, header=True)

for chr in CHROMOSOMES:
	df_chrom = df[df['#CHROM'] == chr]
	df_chrom = df_chrom.sort_values(by='POS')
	df_chrom[['#CHROM','POS','ID','REF','ALT']].to_csv(args.output_file_prefix + '_chr' + chr + '.vcf', sep='\t', index=False, header=True)
