from modules.common import read_fasta
import argparse

parser = argparse.ArgumentParser(description='Reads a FASTA file, removes entries only matching to UTR regions.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input FASTA file", metavar="FILE")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output FASTA file", metavar="FILE")

args = parser.parse_args()

all_proteins = read_fasta(args.input_file)

outfile = open(args.output_file, 'w')

for protein in all_proteins.values():
	proteinIDs = protein['description'].split('protein_IDs:', 1)[1].split(maxsplit=1)[0].split(';')
	is_UTR = all([ 'UTR' in id for id in proteinIDs ])
	if not is_UTR:
		outfile.write('>' + protein['tag'] + '|' + protein['accession'] + '|' + protein['description'] + '\n')
		outfile.write(protein['sequence'] + '\n')

outfile.close()
