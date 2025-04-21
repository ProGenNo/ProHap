from modules.common import read_fasta, KeyWrapper
import bisect
import argparse
import gzip

parser = argparse.ArgumentParser(description='Reads a FASTA file, aggregated all the duplicate sequences into one entry.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input FASTA file", metavar="FILE")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output FASTA file", metavar="FILE")

args = parser.parse_args()

all_proteins = read_fasta(args.input_file)

result_proteins = []

for protein in all_proteins.values():
	seq = protein['sequence']
	seq_hash = hash(seq)

	if ('hap' in protein['tag']) or ('var' in protein['tag']):
		matching_proteins = protein['description'].split('matching_proteins:', 1)[1].split(maxsplit=1)[0].split(';')
		rfs = protein['description'].split('reading_frame:', 1)[1].split(maxsplit=1)[0].split(';')
		protein_start = protein['description'].split('start:', 1)[1].split(maxsplit=1)[0]
		if ('position_within_protein' in protein['description']):
			seq_position = protein['description'].split('position_within_protein:', 1)[1].split(maxsplit=1)[0]
		else:
			seq_position = '0'
	elif 'ref' in protein['tag']:
		matching_proteins = protein['description'].split('matching_proteins:', 1)[1].split(maxsplit=1)[0].split(';')
		rfs = ['-']
		protein_start = '0'
		if ('position_within_protein' in protein['description']):
			seq_position = protein['description'].split('position_within_protein:', 1)[1].split(maxsplit=1)[0]
		else:
			seq_position = '0'
	else:
		matching_proteins = [protein['accession']]
		rfs = ['-']
		protein_start = '0'
		seq_position = '0'

	nearest_idx = bisect.bisect_left(KeyWrapper(result_proteins, key=lambda x: x['hash']), seq_hash)
	if (len(result_proteins) > nearest_idx and result_proteins[nearest_idx]['hash'] == seq_hash):
		result_proteins[nearest_idx]['matching_proteins'].append(matching_proteins)
		result_proteins[nearest_idx]['rfs'].append(rfs)
		result_proteins[nearest_idx]['split_sequences'].append(protein['accession'])
		result_proteins[nearest_idx]['tags'].append(protein['tag'])
		result_proteins[nearest_idx]['start'].append(protein_start)
		result_proteins[nearest_idx]['seq_position'].append(seq_position)
	else:
		result_proteins.insert(nearest_idx, {'hash': seq_hash, 'matching_proteins': [matching_proteins], 'split_sequences': [protein['accession']], 'tags': [protein['tag']], 'sequence': seq, 'start': [protein_start], 'seq_position': [seq_position], 'rfs': [rfs]})

outfile = gzip.open(args.output_file, 'wt') if args.output_file.endswith('.gz') else open(args.output_file, 'w')

for i,protein in enumerate(result_proteins):
	matching_proteins = [ ','.join(plist) for plist in protein['matching_proteins'] ]
	rfs = [ ','.join(rflist) for rflist in protein['rfs'] ]
	description = 'position_within_protein:' + ';'.join(protein['seq_position']) + ' start:' + ';'.join(protein['start']) + ' matching_proteins:' + ';'.join(matching_proteins) + ' reading_frame:' + ';'.join(rfs) + ' split_sequences:' + ';'.join(protein['split_sequences'])

	tag = 'generic_'
	if ('generic_cont' in protein['tags']) or ('generic_sp' in protein['tags']):
		tag += 'cont'
	elif ('generic_ensref' in protein['tags']):
		tag += 'ensref'
	elif ('generic_ensvar' in protein['tags']):
		tag += 'ensvar'
	elif ('generic_var' in protein['tags']):
		tag += 'var'
	elif ('generic_manual' in protein['tags']):		# in case some sequences were added to the fasta manually
		tag += 'manual'
	elif ('generic_enshap' in protein['tags']):
		tag += 'enshap'
	elif ('generic_decoyvar' in protein['tags']):
		tag += 'decoyvar'
	else:
		tag += 'other'

	outfile.write('>' + tag + '|prot_' + hex(i)[2:] + '|' + description.replace('\n', '') + '\n')
	outfile.write(protein['sequence'] + '\n')

outfile.close()
