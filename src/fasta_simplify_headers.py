from modules.common import read_fasta
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Reads a FASTA file, extracts all extra information in headers into a separate tab-separated file.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input FASTA file", metavar="FILE")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output FASTA file", metavar="FILE")

parser.add_argument("-h", dest="output_file_header", required=True,
                    help="output header file (tab-separated)", metavar="FILE")

args = parser.parse_args()

all_proteins = read_fasta(args.input_file)

outfile = open(args.output_file, 'w')
header_data = []

for protein in all_proteins.values():
    # extract all necessary info from the header
    proteinIDs = protein['description'].split('matching_proteins:', 1)[1].split(maxsplit=1)[0]
    proteinPos = protein['description'].split('position_within_protein:', 1)[1].split(maxsplit=1)[0]
    proteinStart = protein['description'].split('start:', 1)[1].split(maxsplit=1)[0]
    proteinRF = protein['description'].split('reading_frame:', 1)[1].split(maxsplit=1)[0]
	
    header_data.append([protein['accession'], protein['tag'], proteinIDs, proteinPos, proteinStart, proteinRF])

    # write the simplified fasta entry
    outfile.write('>' + protein['accession'] + '\n')
    outfile.write(protein['sequence'] + '\n')

outfile.close()

header_df = pd.DataFrame(data=header_data, columns=['accession', 'tag', 'matching_proteins', 'position_within_protein', 'start', 'reading_frame'])
header_df.to_csv(args.output_file_header, sep='\t', index=False)