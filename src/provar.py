import gffutils
import argparse
import os
import pandas as pd
from datetime import datetime

from modules.vcf_reader import parse_vcf
from modules.common import check_open_file, read_fasta
from modules.process_variants import process_store_variants, empty_output

parser = argparse.ArgumentParser(
        description='Creates a database and a FASTA file of variant protein sequences given a VCF file')

parser.add_argument("-i", dest="input_vcf", required=True,
                    help="input VCF file", metavar="FILE",
                    type=lambda x: check_open_file(parser, x))

parser.add_argument("-db", dest="annotation_db", required=True,
                    help="DB file created by gffutils from GTF")

parser.add_argument("-af", dest="min_af", required=False, type=float,
                    help="Allele Frequency (AF) lower threshold - default 0", default=0)

parser.add_argument("-cdna", dest="cdnas_fasta", required=True,
                    help="input cDNA fasta file")

parser.add_argument("-chr", dest="chromosome", required=True,
                    help="chromosome being processed (e.g., 1, 12 or X)")

parser.add_argument("-transcripts", dest="transcript_list", required=True,
                    help="list of transcript IDs, provided in a CSV file", metavar="FILE")

parser.add_argument("-require_start", dest="require_start", required=False, type=int,
                    help="flag: require annotation of the start codon, set to 0 to disable; default: 1", default=1)

parser.add_argument("-force_rf", dest="force_rf", required=False,
                    help="Force the most likely reading frame when start codon is not annotated or lost due to mutation, set to 0 to disable; default: 1", default=1)

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

parser.add_argument("-output_cdna_fasta", dest="output_cdna_fasta", required=False, default="",
                    help="output cDNA FASTA file (optional; default: none)")

args = parser.parse_args()

print('[ProVar] Generating protein variants from', args.input_vcf.name)

print (('Chr ' + args.chromosome + ':'), 'Reading', args.annotation_db)
# Load the annotations database
annotations_db = gffutils.FeatureDB(args.annotation_db)

print (('Chr ' + args.chromosome + ':'), 'Reading', args.transcript_list)
# read the list of transcript IDs
transcript_df = pd.read_csv(args.transcript_list)
transcript_df['chromosome'] = transcript_df['chromosome'].apply(lambda x: str(x))
transcript_list = transcript_df[transcript_df['chromosome'] == args.chromosome]['transcriptID'].tolist()

print (('Chr ' + args.chromosome + ':'), 'Assigning annotations to transcripts.')

# create a list of transcript features from the annotation db
all_transcripts = []

for transcript_id in transcript_list:
    feature = annotations_db[transcript_id]
    if (args.require_start):
        start_codons = [ sc for sc in annotations_db.children(feature, featuretype='start_codon', order_by='start') ]    # there should be only one, but just in case...
        if (len(start_codons) > 0):
                all_transcripts.append(feature)
    else:
        all_transcripts.append(feature)

all_transcripts.sort(key=lambda x: x.start)
transcript_list = [ feature.id for feature in all_transcripts ]

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
        all_cds = read_fasta(args.cdnas_fasta, True)

        log_file = open(args.log_file, 'a')
        log_file.write('------------' + '[' + datetime.now().strftime('%X %x') + '] Chr ' + args.chromosome + ':' + '------------\n')

        print (('Chr ' + args.chromosome + ':'), 'Creating variant database.')
        # align the variant coordinates to transcript, translate into the protein database
        process_store_variants(all_transcripts, args.tmp_dir, log_file, all_cds, annotations_db, args.chromosome, args.fasta_tag, args.accession_prefix, args.force_rf, args.output_file, args.output_fasta, args.output_cdna_fasta)

        log_file.close()

        # remove the temporary files
        for transcript_id in transcript_list:
                os.remove(args.tmp_dir + '/' + transcript_id + '.tsv')

        print('Done.')
