'''
Main file of the ProHap tool

Creates a database of CDS + protein haplotypes, and a fasta file of protein haplotype sequences.
'''

import gffutils
import argparse
import os
from numpy import int64
import pandas as pd

from modules.vcf_reader import parse_vcf
from modules.common import read_fasta
from modules.get_haplotypes import get_gene_haplotypes
from modules.process_haplotypes import process_store_haplotypes, empty_output

parser = argparse.ArgumentParser(
        description='Creates a database of CDS + protein haplotypes, and a fasta file of protein haplotype sequences.')

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

parser.add_argument("-s", dest="samples_filename", required=True,
                    help="tab-separated file with sample information (must include 'Sample name' and 'Sex' columns)")

parser.add_argument("-af", dest="min_af", required=False, type=float,
                    help="Allele Frequency (AF) lower threshold - default 0", default=0)

parser.add_argument("-cdna", dest="cdnas_fasta", required=True,
                    help="input cDNA fasta file")

parser.add_argument("-transcripts", dest="transcript_list", required=True,
                    help="list of transcript IDs, provided in a file", metavar="FILE",
                    type=lambda x: is_valid_file(parser, x))

parser.add_argument("-require_start", dest="require_start", required=False, type=int,
                    help="flag: require annotation of the start codon, set to 0 to disable; default: 1", default=1)

parser.add_argument("-x_par1_to", dest="x_par1_to", required=False, type=int64,
                    help="end location of the 1st pseudoautosomal region on chromosome X; default: 2,781,479", default=2781479)

parser.add_argument("-x_par2_from", dest="x_par2_from", required=False, type=int64,
                    help="start location of the 2nd pseudoautosomal region on chromosome X; default: 155,701,383", default=155701383)

parser.add_argument("-chr", dest="chromosome", required=True,
                    help="chromosome being processed (e.g., 1, 12 or X)")

parser.add_argument("-threads", dest="threads", required=False, type=int,
                    help="number of threads to use; default: 4", default=4)

parser.add_argument("-foo", dest="min_foo", required=False, type=float,
                    help="Minimum FoO of a haplotype to be reported in the fasta", default=0)
                    
parser.add_argument("-tag", dest="fasta_tag", required=False,
                    help="tag for FASTa file entries", default='generic_enshap')

parser.add_argument("-id_prefix", dest="haplo_id_prefix", required=False,
                    help="prefix for the haplotype identifier", default='haplo_')

parser.add_argument("-acc_prefix", dest="accession_prefix", required=False,
                    help="prefix for FASTA file entries accession", default='enshap')

parser.add_argument("-log", dest="log_file", required=False,
                    help="output log file", default="prohap.log")

parser.add_argument("-tmp_dir", dest="tmp_dir", required=False,
                    help="directory for temporary files", default="tmp")

parser.add_argument("-output_csv", dest="output_file", required=True,
                    help="output CSV file")

parser.add_argument("-output_fasta", dest="output_fasta", required=True,
                    help="output FASTA file")

args = parser.parse_args()

print('[ProHap] Computing protein haplotypes from', args.input_vcf.name)

print (('Chr ' + args.chromosome + ':'), 'Reading', args.annotation_db)
# Load the annotations database
annotations_db = gffutils.FeatureDB(args.annotation_db)

print (('Chr ' + args.chromosome + ':'), 'Reading', args.transcript_list.name)
# read the list of transcript IDs
transcript_list = [ line[:-1].split('.', 1)[0] for line in args.transcript_list.readlines() ]

print (('Chr ' + args.chromosome + ':'), 'Reading', args.samples_filename)
# get sample IDs of males from the metadata file
samples_df = pd.read_csv(args.samples_filename, sep='\t')
male_samples = samples_df[samples_df['Sex'] == 'male']['Sample name'].tolist()

print (('Chr ' + args.chromosome + ':'), 'Assigning annotations to transcripts.')
# create a list of transcript features from the annotation db
all_transcripts = []
for transcript_id in transcript_list:
    all_transcripts.append(annotations_db[transcript_id])

# check if start codon is annotated
if (args.require_start):
        filtered_features = []

        for feature in all_transcripts:
                start_codons = [ sc for sc in annotations_db.children(feature, featuretype='start_codon', order_by='start') ]    # there should be only one, but just in case...
                if (len(start_codons) > 0):
                        filtered_features.append(feature)

        all_transcripts = filtered_features
        transcript_list = [ feature.id for feature in all_transcripts ]

all_transcripts.sort(key=lambda x: x.start)

print (('Chr ' + args.chromosome + ':'), 'Assigning variants to transcripts.')
# parse the VCF file, get a dataframe of variants for each transcript
vcf_colnames = parse_vcf(all_transcripts, args.input_vcf, annotations_db, args.min_af, args.tmp_dir)
#vcf_colnames = list(pd.read_csv(args.tmp_dir + '/' + transcript_list[0] + '.tsv', sep='\t').columns.values)

# check if the vcf file was empty
if (len(vcf_colnames) == 0):
        print(('Chr ' + args.chromosome + ':'), 'VCF file is empty, creating empty output files.')
        empty_output(args.output_file, args.output_fasta)
else:
        print (('Chr ' + args.chromosome + ':'), 'Computing co-occurence of alleles.')
        # check co-occurence of alleles -> get the haplotypes for all transcripts
        gene_haplo_df = get_gene_haplotypes(all_transcripts, vcf_colnames, args.tmp_dir, args.log_file, args.threads, (args.chromosome == 'X'), args.x_par1_to, args.x_par2_from, male_samples)

        # remove the temporary files
        for transcript_id in transcript_list:
                os.remove(args.tmp_dir + '/' + transcript_id + '.tsv')

        # filter the haplotypes by FoO
        gene_haplo_df = gene_haplo_df[gene_haplo_df['Frequency'] >= args.min_foo]

        # read the CDS sequence file
        print (('Chr ' + args.chromosome + ':'), "Reading", args.cdnas_fasta)
        all_cds = read_fasta(args.cdnas_fasta)

        print (('Chr ' + args.chromosome + ':'), 'Creating haplotype database.')
        # align the variant coordinates to transcript, translate into the protein database
        process_store_haplotypes(gene_haplo_df, all_cds, annotations_db, args.chromosome, args.fasta_tag, args.haplo_id_prefix, args.accession_prefix, args.output_file, args.output_fasta)
