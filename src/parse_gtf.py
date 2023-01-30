import gffutils
import argparse

parser = argparse.ArgumentParser(
        description='Reads the GTF file, parses it into the sqlite3-like DB file used by gffutils.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input GTF")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output DB")

args = parser.parse_args()

# parse the GTF file in to the database -> write the DB file
db = gffutils.create_db(args.input_file, dbfn=args.output_file, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True, disable_infer_genes=True, disable_infer_transcripts=True)
