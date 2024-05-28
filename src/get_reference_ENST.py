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
    chr = '-'
    if ('chromosome:' in prot['description']):
        chr = prot['description'].split('chromosome:')[1].split(':',2)[1]
    elif ('GRCh' in prot['description']):
        chr = prot['description'].split('GRCh',1)[1].split(':',2)[1]  

    result_data.append([chr,trID])

result_df = pd.DataFrame(data=result_data, columns=['chromosome', 'transcriptID'])

if (args.only_MANE):    
        # check for which of the included genes there is a MANE Select transcript available
        result_df['geneID'] = result_df['transcriptID'].apply(lambda trID: annotations_db[trID].attributes['gene_id'][0])
        gene_data = [ [ geneID, any([(('tag' in tr.attributes) and ('MANE_Select' in tr.attributes['tag'])) for tr in annotations_db.children(annotations_db[geneID], featuretype='transcript')]) ] for geneID in result_df['geneID'].drop_duplicates().tolist() ]
        genes_df = pd.DataFrame(data=gene_data, columns=['geneID', 'has_MANE']).set_index('geneID')

        transcript_filter = []      # mask to filter out non-canonical transcripts - contains a boolean value for each row of result_df
        
        for i,row in result_df.iterrows():
                tr_feature = annotations_db[row['transcriptID']]

                # if there is a MANE select transcript for this gene, keep only that transcript, otherwise keep the Ensembl Canonical transcript
                if (genes_df.loc[row['geneID']]['has_MANE']):
                       transcript_filter.append(('tag' in tr_feature.attributes) and ('MANE_Select' in tr_feature.attributes['tag']))
                else:
                       transcript_filter.append(('tag' in tr_feature.attributes) and ('Ensembl_canonical' in tr_feature.attributes['tag'])) 

        result_df.drop('geneID', axis=1, inplace=True)
        result_df = result_df[transcript_filter]

result_df.sort_values(by='chromosome', inplace=True)
result_df.to_csv(args.output_file, index=False)