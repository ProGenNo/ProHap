import gffutils
import argparse
import os
from datetime import datetime

from vcf_reader import parse_vcf
from common import read_fasta
from process_variants import process_store_variants, empty_output

parser = argparse.ArgumentParser(
        description='Creates a database and a FASTA file of variant protein sequences given a VCF file')

def is_valid_file(parser, arg):
        if not os.path.exists(arg):
                parser.error("The file %s does not exist!" % arg)
        else:
                return open(arg, 'r')  # return an open file handle

parser.add_argument("-i", dest="input_vcf", required=True,
                    help="input VCF file", metavar="FILE",
                    type=lambda x: is_valid_file(parser, x))

parser.add_argument("-db", dest="annotation_db", required=True,
                    help="DB file created by gffutils from GTF")

parser.add_argument("-af", dest="min_af", required=False, type=float,
                    help="Allele Frequency (AF) lower threshold - default 0", default=0)

parser.add_argument("-cdna", dest="cdnas_fasta", required=True,
                    help="input cDNA fasta file")

parser.add_argument("-chr", dest="chromosome", required=True,
                    help="chromosome being processed (e.g., 1, 12 or X)")

parser.add_argument("-transcripts", dest="transcript_list", required=True,
                    help="list of transcript IDs, provided in a file", metavar="FILE",
                    type=lambda x: is_valid_file(parser, x))

parser.add_argument("-tag", dest="fasta_tag", required=False,
                    help="tag for FASTa file entries", default='generic_var')

parser.add_argument("-acc_prefix", dest="accession_prefix", required=False,
                    help="prefix for FASTA file entries accession", default='var')

parser.add_argument("-log", dest="log_file", required=False,
                    help="output log file", default="provar.log")

parser.add_argument("-tmp_dir", dest="tmp_dir", required=False,
                    help="directory for temporary files", default="tmp")

parser.add_argument("-output_csv", dest="output_file", required=True,
                    help="output CSV file")

parser.add_argument("-output_fasta", dest="output_fasta", required=True,
                    help="output FASTA file")

args = parser.parse_args()

print('[ProVar] Generating protein variants from', args.input_vcf.name)

print (('Chr ' + args.chromosome + ':'), 'Reading', args.annotation_db)
# Load the annotations database
annotations_db = gffutils.FeatureDB(args.annotation_db)

print (('Chr ' + args.chromosome + ':'), 'Reading', args.transcript_list.name)
# read the list of transcript IDs
transcript_list = [ line[:-1].split('.', 1)[0] for line in args.transcript_list.readlines() ]

print (('Chr ' + args.chromosome + ':'), 'Assigning annotations to transcripts.')
# create a list of transcript features from the annotation db
all_transcripts = []
for transcript_id in transcript_list:
    all_transcripts.append(annotations_db[transcript_id])

all_transcripts.sort(key=lambda x: x.start)

print (('Chr ' + args.chromosome + ':'), 'Assigning variants to transcripts.')
# parse the VCF file, get a dataframe of variants for each transcript
vcf_columns = parse_vcf(all_transcripts, args.input_vcf, annotations_db, args.min_af, args.tmp_dir)

# check if the vcf file was empty
if (len(vcf_columns) == 0):
        print(('Chr ' + args.chromosome + ':'), 'VCF file is empty, creating empty output files.')
        empty_output(args.output_file, args.output_fasta)
else:
        # read the CDS sequence file
        print (('Chr ' + args.chromosome + ':'), "Reading", args.cdnas_fasta)
        all_cds = read_fasta(args.cdnas_fasta)

        log_file = open(args.log_file, 'a')
        log_file.write('------------' + '[' + datetime.now().strftime('%X %x') + '] Chr ' + args.chromosome + ':' + '------------\n')

        print (('Chr ' + args.chromosome + ':'), 'Creating variant database.')
        # align the variant coordinates to transcript, translate into the protein database
        process_store_variants(all_transcripts, args.tmp_dir, log_file, all_cds, annotations_db, args.chromosome, args.fasta_tag, args.accession_prefix, args.output_file, args.output_fasta)

        log_file.close()

        # remove the temporary files
        for transcript_id in transcript_list:
                os.remove(args.tmp_dir + '/' + transcript_id + '.tsv')

        print('Done.')