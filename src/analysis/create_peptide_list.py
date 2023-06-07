import json
import pandas as pd
import argparse
import bisect
from common import read_fasta, digest

parser = argparse.ArgumentParser(description="Creates a list of discoverable peptides from a FASTA database")

parser.add_argument("-i", dest="input_file", required=True,
                    help="input FASTA file", metavar="FILE")

parser.add_argument("-m", dest="missed_cl", type=int, required=False, default=2,
                    help="# missed cleavage sites allowed (default: 2")

parser.add_argument("-enz", dest="enzyme", required=False, default="Trypsin",
                    help="Name of the simulated proteases cleavage pattern (default: \"Trypsin\"")

parser.add_argument("-enz_json", dest="enzyme_json", required=False, default="data/enzymeFactory-5.0.0.json",
                    help="JSON file with enzyme information", metavar="FILE")

parser.add_argument("-min_len", dest="min_len", type=int, required=False, default=6,
                    help="minimum peptide length (default: 6")

parser.add_argument("-max_len", dest="max_len", type=int, required=False, default=40,
                    help="maximum peptide length (default: 40")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output TSV file", metavar="FILE")

args = parser.parse_args()

print("Reading file", args.input_file)
all_proteins = read_fasta(args.input_file)

json_file = open(args.enzyme_json, 'r')
enzyme_patterns = json.load(json_file)
if (args.enzyme not in enzyme_patterns['enzymes']):
    raise Exception("Cleavage pattern for the enzyme " + args.enzyme + " not found in the provided JSON: " + args.enzyme_json)
cleavage_pattern = enzyme_patterns['enzymes'][args.enzyme]
json_file.close()

all_peptides = []
all_matching_proteins = []
all_peptide_positions = []

print ('Creating the peptide list...')
for prot_idx,protein in enumerate(list(all_proteins.values())):
    peptides, peptide_posotions = digest(protein['sequence'], args.missed_cl, args.min_len, args.max_len, cleavage_pattern)

    for i,peptide in enumerate(peptides):
        idx = bisect.bisect_left(all_peptides, peptide)

        if ((idx >= len(all_peptides)) or (all_peptides[idx] != peptide)):
            all_peptides.insert(idx, peptide)
            all_matching_proteins.insert(idx, [protein['accession']])
            all_peptide_positions.insert(idx, [str(peptide_posotions[i])])
        else:
            all_matching_proteins[idx].append(protein['accession'])
            all_peptide_positions[idx].append(str(peptide_posotions[i]))
    print (prot_idx, '/', len(all_proteins), end='\r')


df_data = [ ['pep_' + hex(i)[2:], all_peptides[i], ';'.join(all_matching_proteins[i]), ';'.join(positions)] for i,positions in enumerate(all_peptide_positions) ]

output_df = pd.DataFrame(data=df_data, columns=['ID', 'Sequence', 'Proteins', 'Positions'])
output_df.to_csv(args.output_file, header=True, index=False, sep='\t')
