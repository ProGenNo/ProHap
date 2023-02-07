import argparse
import pandas as pd
import re
import vcf

parser = argparse.ArgumentParser(description="Splits a VCF into separate files for each chromosome, fills in missing values by '-'.")

parser.add_argument("-i", dest="input_file", required=True,
                    help="input CSV file", metavar="FILE")

parser.add_argument("-o", dest="output_file_prefix", required=True,
                    help="output VCF file prefix", metavar="FILE")

args = parser.parse_args()

# if (not re.match(r'[CGTA]', str(ALT)))

infile = open(args.input_file, 'r')
vcf_reader = vcf.Reader(infile)
vcf_data = []

def fill_missing(val):
	if (not re.match(r'[CGTA]', str(val))):
		return '-'
	return val

for record in vcf_reader:
    for i,alt in enumerate(record.ALT):        
        MAF = -1
        if ('AF' in record.INFO):
            MAF = record.INFO['AF'][i]
        elif ('MAF' in record.INFO):
            MAF = record.INFO['MAF'][i]

        ID = '.'
        if record.ID is not None:
            ID = record.ID

        vcf_data.append([record.CHROM.replace('chr', ''), record.POS, ID, fill_missing(record.REF), fill_missing(alt.sequence), 'MAF=' + str(MAF)])

df = pd.DataFrame(data=vcf_data, columns=['#CHROM','POS','ID','REF','ALT','INFO'])
CHROMOSOMES = [str(x) for x in list(range(1, 23))] + ['X', 'Y']

for chr in CHROMOSOMES:
	df_chrom = df[df['#CHROM'] == chr]
	df_chrom = df_chrom.sort_values(by='POS')
	df_chrom.to_csv(args.output_file_prefix + '_chr' + chr + '.vcf', sep='\t', index=False, header=True)

infile.close()