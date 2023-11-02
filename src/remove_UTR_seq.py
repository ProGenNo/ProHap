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
	splitIDs = protein['description'].split('split_sequences:', 1)[1].split(maxsplit=1)[0].split(';')

	# remove information for proteins that contain only their UTR region in this sequence
	proteinIDs = [ x for i,x in enumerate(protein['description'].split('matching_proteins:', 1)[1].split(maxsplit=1)[0].split(';')) if ('UTR' not in splitIDs[i]) ]
	proteinPos = [ x for i,x in enumerate(protein['description'].split('position_within_protein:', 1)[1].split(maxsplit=1)[0].split(';')) if ('UTR' not in splitIDs[i]) ]
	proteinStart = [ x for i,x in enumerate(protein['description'].split('start:', 1)[1].split(maxsplit=1)[0].split(';')) if ('UTR' not in splitIDs[i]) ]
	proteinRF = [ x for i,x in enumerate(protein['description'].split('reading_frame:', 1)[1].split(maxsplit=1)[0].split(';')) if ('UTR' not in splitIDs[i]) ]

	if (len(proteinIDs) > 0):
		# remove the split sequence IDs from the description - they are redundant for the search and PSM annotation
		description = 'position_within_protein:' + ';'.join(proteinPos) + ' start:' + ';'.join(proteinStart) + ' matching_proteins:' + ';'.join(proteinIDs) + ' reading_frame:' + ';'.join(proteinRF)
		outfile.write('>' + protein['tag'] + '|' + protein['accession'] + '|' + description + '\n')
		outfile.write(protein['sequence'] + '\n')

outfile.close()
