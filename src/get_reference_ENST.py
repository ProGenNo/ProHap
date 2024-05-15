import argparse
import pandas as pd
import gffutils
from modules.common import read_fasta

parser = argparse.ArgumentParser(
        description='Reads the reference FASTA file, outputs a list of stable IDs for transcripts that are included in this database.')

parser.add_argument("-i", dest="ref_fasta", required=True,
                    help="input fasta file")

parser.add_argument("-annot", dest="annotation_db", required=True,
                    help="annotations DB file")

parser.add_argument("-MANE", dest="only_MANE", required=False, type=int,
                    help="include only MANE Select transcripts (default: 0)", default=0)

parser.add_argument("-o", dest="output_file", required=True,
                    help="output transcript list file")

args = parser.parse_args()

print ("Reading", args.ref_fasta)
ref_proteins = read_fasta(args.ref_fasta)

# Load the annotations database
print ('Reading', args.annotation_db)
annotations_db = gffutils.FeatureDB(args.annotation_db)

CHROMOSOMES = [str(x) for x in list(range(1, 23))] + ['X', 'Y']

result_data = []

for prot in ref_proteins.values():
    trID = prot['description'].split('transcript:',1)[1].split('.',1)[0]
    chr = prot['description'].split(':',3)[2]

    if chr in CHROMOSOMES:
        transcript_feature = annotations_db[trID]
        if (('tag' in transcript_feature.attributes) and ('MANE_Select' in transcript_feature.attributes['tag'])):
                result_data.append([chr,trID])

result_df = pd.DataFrame(data=result_data, columns=['chromosome', 'transcriptID'])
result_df.sort_values(by='chromosome', inplace=True)
result_df.to_csv(args.output_file, index=False)