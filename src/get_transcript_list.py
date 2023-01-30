import gffutils
import argparse

parser = argparse.ArgumentParser(
        description='Reads the DB file created by gffutils, outputs a list of stable IDs for transcripts that match listed biotype annotations.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input DB file")

parser.add_argument('-bio', 'biotypes', required=False, type=str,
                    help="List of permitted transcript biotypes, comma-separated (default: protein_coding). Specify \"all\" to disable filtering.", default="protein_coding")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output transcript list file")

args = parser.parse_args()

# check if at least one of the annotated biotypes is in the list
def check_transcript(transcript_feature, valid_biotypes):
    filter_pass = False
    for biotype in transcript_feature['transcript_biotype']:
        filter_pass = filter_pass or (biotype in valid_biotypes)

    return filter_pass

# Load the annotations database
print ('Reading', args.annotation_db)
annotations_db = gffutils.FeatureDB(args.annotation_db)

biotypes_filter = args.biotypes.split(',')

print ('Filtering transcripts')
# loop through all the transcript features and filter
transcript_features = [ feature for feature in annotations_db.features_of_type('transcript', order_by='start') if ((args.biotypes == 'all') or check_transcript(feature, biotypes_filter))]

print ('Writing', args.output_file)
# write into the file
outfile = open(args.output_file, 'w')
outfile.write("chromosome,transcriptID\n")  # header
for transcript in transcript_features:
        outfile.write(transcript.chrom + ',' + transcript.id + '\n')
outfile.close()