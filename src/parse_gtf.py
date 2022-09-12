import gffutils
import argparse

parser = argparse.ArgumentParser(
        description='Reads the GTF file, parses it into the sqlite3-like DB file used by gffutils. Outputs a list of stable IDs for transcripts that satisfy given criteria.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input GTF")

parser.add_argument("-noncoding", dest="include_noncoding", required=False, type=int,
                    help="flag: include non-coding transcripts; default: 0", default=0)

parser.add_argument("-o", dest="output_file", required=True,
                    help="output DB")

parser.add_argument("-transcript_list", dest="transcript_list", required=True,
                    help="output transcript list file")

args = parser.parse_args()

# parse the GTF file in to the database -> write the DB file
db = gffutils.create_db(args.input_file, dbfn=args.output_file, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True, disable_infer_genes=True, disable_infer_transcripts=True)

# extraxt transcript IDs, filter only protein-coding if necessary
transcript_features = [ feature for feature in db.features_of_type('transcript', order_by='start') if (args.include_noncoding or ('protein_coding' in feature['transcript_biotype']))]

# write into the file
outfile = open(args.transcript_list, 'w')
for transcript in transcript_features:
        outfile.write(transcript.id + '\n')
outfile.close()