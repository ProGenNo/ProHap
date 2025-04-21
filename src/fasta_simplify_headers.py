from modules.common import read_fasta
import pandas as pd
import gzip
import re
import argparse
import gffutils

parser = argparse.ArgumentParser(description='Reads a FASTA file, extracts all extra information in headers into a separate tab-separated file.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input FASTA file", metavar="FILE")

parser.add_argument("-hap_tsv", dest="haplo_table", required=False,
                    help="input haplotype table", metavar="FILE")

parser.add_argument("-hap_prefix", dest="haplo_prefix", required=False,
                    help="prefix for haplotype protein ID (default: 'haplo_')", default='haplo_')

parser.add_argument("-var_prefix", dest="var_prefix", required=False,
                    help="prefix for variant protein ID (default: 'var_')", default='var_')

parser.add_argument("-db", dest="annotation_db", required=True,
                    help="DB file created by gffutils from GTF")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output FASTA file", metavar="FILE")

parser.add_argument("-header", dest="output_file_header", required=True,
                    help="output header file (tab-separated)", metavar="FILE")

args = parser.parse_args()

all_proteins = read_fasta(args.input_file)

haplo_df = None
if (args.haplo_table):
    haplo_df = pd.read_table(args.haplo_table)
    haplo_df.set_index('HaplotypeID', inplace=True)

annot = gffutils.FeatureDB(args.annotation_db)

outfile = gzip.open(args.output_file, 'wt') if args.output_file.endswith('.gz') else open(args.output_file, 'w')
header_data = []

for protein in all_proteins.values():
    # extract all necessary info from the header
    proteinIDs = protein['description'].split('matching_proteins:', 1)[1].split(maxsplit=1)[0]
    proteinPos = protein['description'].split('position_within_protein:', 1)[1].split(maxsplit=1)[0]
    proteinStart = protein['description'].split('start:', 1)[1].split(maxsplit=1)[0]
    proteinRF = protein['description'].split('reading_frame:', 1)[1].split(maxsplit=1)[0]

    # if the protein is a contaminant, mark it and don't look for a gene name (use the UniProt contaminant name)
    if (protein['tag'] == 'generic_cont'):
        contIDs = [ protID for protID in re.split(r'[;,]', proteinIDs) if not (protID.startswith('ENST') or protID.startswith(args.haplo_prefix) or protID.startswith(args.var_prefix)) ]        
        outfile.write('>' + protein['accession'] + ' CONTAMINANT GN=' + ';'.join(contIDs) + '\n')
        outfile.write(protein['sequence'] + '\n')
        header_data.append([protein['accession'], protein['tag'], proteinIDs, proteinPos, proteinStart, proteinRF])
        continue

    # aggregate all the matching transcripts for this protein
    if (protein['tag'] == 'generic_ensref'):
        trIDs = [ protID for protID in re.split(r'[;,]', proteinIDs) if (protID.startswith('ENST')) ]
    elif (protein['tag'] == 'generic_enshap'):
        trIDs = [ haplo_df.loc[protID]['TranscriptID'] for protID in re.split(r'[;,]', proteinIDs) if (protID.startswith(args.haplo_prefix)) ]
    elif (protein['tag'] == 'generic_var'):
        trIDs = [ 'ENST' + protID.split('ENST',1)[1].split('_',1)[0] for protID in re.split(r'[;,]', proteinIDs) if (protID.startswith(args.var_prefix)) ]
    else:
        trIDs = []
    
    # based on the list of transcripts, get the corresponding gene name (or the ENSG identifier if name not available)
    gIDs = []
    for trID in trIDs:
        feature_attr = annot[trID].attributes
        if ('gene_id' in feature_attr):
            gIDs.extend(feature_attr['gene_name'] if ('gene_name' in feature_attr) else feature_attr['gene_id'])

    gIDs = list(dict.fromkeys(gIDs))

    # write the simplified fasta entry
    if (len(gIDs) > 1):
        # if the protein matches more than one gene, write a duplicate entry for each to avoid confusion
        for i,gID in enumerate(gIDs):
            outfile.write('>' + protein['accession'] + '.' + str(i) + ' GN=' + gID + '\n')
            outfile.write(protein['sequence'] + '\n')	
            header_data.append([protein['accession'] + '.' + str(i), protein['tag'], proteinIDs, proteinPos, proteinStart, proteinRF])

    elif (len(gIDs) == 1):
        outfile.write('>' + protein['accession'] + ' GN=' + gIDs[0] + '\n')
        outfile.write(protein['sequence'] + '\n')
        header_data.append([protein['accession'], protein['tag'], proteinIDs, proteinPos, proteinStart, proteinRF])
        
    else:
        outfile.write('>' + protein['accession'] + ' GN=missing\n')
        outfile.write(protein['sequence'] + '\n')
        header_data.append([protein['accession'], protein['tag'], proteinIDs, proteinPos, proteinStart, proteinRF])

outfile.close()

header_df = pd.DataFrame(data=header_data, columns=['accession', 'tag', 'matching_proteins', 'position_within_protein', 'start', 'reading_frame'])
header_df.to_csv(args.output_file_header, sep='\t', index=False, compression='infer')