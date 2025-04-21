'''
Main file of the ProHap tool

Creates a database of CDS + protein haplotypes, and a fasta file of protein haplotype sequences.
'''

import gffutils
import argparse
import os
import gzip
from numpy import int64
import pandas as pd

from modules.vcf_reader import parse_vcf
from modules.common import check_open_file, read_fasta
from modules.get_haplotypes import get_gene_haplotypes
from modules.process_haplotypes import process_haplotypes, empty_output

parser = argparse.ArgumentParser(
        description='Creates a database of CDS + protein haplotypes, and a fasta file of protein haplotype sequences.')

parser.add_argument("-i", dest="input_vcf", required=True,
                    help="input VCF file", metavar="FILE",
                    type=lambda x: check_open_file(parser, x))

parser.add_argument("-db", dest="annotation_db", required=True,
                    help="DB file created by gffutils from GTF")

parser.add_argument("-s", dest="samples_filename", required=True,
                    help="tab-separated file with sample information (must include 'Sample name' and 'Sex' columns)")

parser.add_argument("-af", dest="min_af", required=False, type=float,
                    help="Allele Frequency (AF) lower threshold - default 0", default=0)

parser.add_argument("-cdna", dest="cdnas_fasta", required=True,
                    help="input cDNA fasta file")

parser.add_argument("-transcripts", dest="transcript_list", required=True,
                    help="list of transcript IDs, provided in a CSV file", metavar="FILE")

parser.add_argument("-require_start", dest="require_start", required=False, type=int,
                    help="flag: require annotation of the start codon, set to 0 to disable; default: 1", default=1)

parser.add_argument("-ignore_UTR", dest="ignore_UTR", required=False, type=int,
                    help="flag (0 or 1): ignore variation in the UTR sequences, do not add UTR translation to proteins; default: 1", default=1)
                    
parser.add_argument("-skip_start_lost", dest="skip_start_lost", required=False, type=int,
                    help="flag (0 or 1): ignore haplotypes where the start codon is lost; default: 1", default=1)

parser.add_argument("-force_rf", dest="force_rf", required=False,
                    help="Force the most likely reading frame when start codon is not annotated or lost due to mutation, set to 0 to disable; default: 1", default=1)

parser.add_argument("-x_par1_to", dest="x_par1_to", required=False, type=int64,
                    help="end location of the 1st pseudoautosomal region on chromosome X; default: 2,781,479", default=2781479)

parser.add_argument("-x_par2_from", dest="x_par2_from", required=False, type=int64,
                    help="start location of the 2nd pseudoautosomal region on chromosome X; default: 155,701,383", default=155701383)

parser.add_argument("-chr", dest="chromosome", required=True,
                    help="chromosome being processed (e.g., 1, 12 or X)")

parser.add_argument("-threads", dest="threads", required=False, type=int,
                    help="number of threads to use; default: 4", default=4)

parser.add_argument("-min_hap_freq", dest="min_freq", required=False, type=float,
                    help="Minimum frequency of a haplotype to be reported in the result (specify -1 to use count threshold instead); default: -1", default=-1)

parser.add_argument("-min_hap_count", dest="min_hap_count", required=False, type=int,
                    help="Minimum count of occurrences of a haplotype to be reported in the result (used only if frequency threshold is not specified); default: 0", default=0)
                    
parser.add_argument("-tag", dest="fasta_tag", required=False,
                    help="tag for FASTa file entries", default='generic_enshap')

parser.add_argument("-id_prefix", dest="haplo_id_prefix", required=False,
                    help="prefix for the haplotype identifier; default: haplo_", default='haplo_')

parser.add_argument("-acc_prefix", dest="accession_prefix", required=False,
                    help="prefix for FASTA file entries accession; default: enshap", default='enshap')

parser.add_argument("-log", dest="log_file", required=False,
                    help="output log file; default: prohap.log", default="prohap.log")

parser.add_argument("-tmp_dir", dest="tmp_dir", required=False,
                    help="directory for temporary files; default: tmp", default="tmp")

parser.add_argument("-output_csv", dest="output_file", required=True,
                    help="output CSV file")

parser.add_argument("-output_fasta", dest="output_fasta", required=True,
                    help="output FASTA file")

parser.add_argument("-output_cdna_fasta", dest="output_cdna_fasta", required=False, default="",
                    help="output cDNA FASTA file (optional; default: none)")

args = parser.parse_args()

print('[ProHap] Computing protein haplotypes from', args.input_vcf.name)

print (('Chr ' + args.chromosome + ':'), 'Reading', args.annotation_db)
# Load the annotations database
annotations_db = gffutils.FeatureDB(args.annotation_db)

print (('Chr ' + args.chromosome + ':'), 'Reading', args.transcript_list)
# read the list of transcript IDs
transcript_df = pd.read_csv(args.transcript_list)
transcript_df['chromosome'] = transcript_df['chromosome'].apply(lambda x: str(x))
transcript_list = transcript_df[transcript_df['chromosome'] == args.chromosome]['transcriptID'].tolist()

print (('Chr ' + args.chromosome + ':'), 'Reading', args.samples_filename)
# get sample IDs of males from the metadata file
samples_df = pd.read_csv(args.samples_filename, sep='\t')
#samples_df.set_index('Sample_name', inplace=True)
#male_samples = samples_df[samples_df['Sex'] == 'male']['Sample name'].tolist()

print (('Chr ' + args.chromosome + ':'), 'Assigning annotations to transcripts.')

# create a list of transcript features from the annotation db
all_transcripts = []

for transcript_id in transcript_list:
    feature = annotations_db[transcript_id]
    if (args.require_start):    # start codon annotation is required - check if present
        start_codons = [ sc for sc in annotations_db.children(feature, featuretype='start_codon', order_by='start') ]    # there should be only one, but just in case...
        if (len(start_codons) > 0):
                all_transcripts.append(feature)
    else:
        all_transcripts.append(feature)

all_transcripts.sort(key=lambda x: x.start)
transcript_list = [ feature.id for feature in all_transcripts ]

print (('Chr ' + args.chromosome + ':'), 'Assigning variants to transcripts.')
# parse the VCF file, get a dataframe of variants for each transcript
vcf_colnames = parse_vcf(all_transcripts, args.input_vcf, annotations_db, args.min_af, args.tmp_dir)

# keep only the samples that are in the metadata table
sample_ids = [sample for sample in vcf_colnames if (sample in samples_df['Sample name'].tolist())]

# check if the vcf file was empty
if (len(vcf_colnames) == 0):
        print(('Chr ' + args.chromosome + ':'), 'VCF file is empty, creating empty output files.')
        empty_output(args.output_file, args.output_fasta)
else:
        print (('Chr ' + args.chromosome + ':'), 'Computing the co-occurrence of alleles.')
        # check co-occurence of alleles -> get the haplotypes for all transcripts
        gene_haplo_df = get_gene_haplotypes(all_transcripts, sample_ids, args.tmp_dir, args.log_file, args.threads, (args.chromosome == 'X'), args.x_par1_to, args.x_par2_from, samples_df)

        # remove the temporary files
        for transcript_id in transcript_list:
                os.remove(args.tmp_dir + '/' + transcript_id + '.tsv')

        # filter the haplotypes by frequency -> CHANGE: filter only after processing, some haplotypes can be merged
        #gene_haplo_df = gene_haplo_df[gene_haplo_df['Frequency'] >= args.min_freq]

        # read the CDS sequence file
        print (('Chr ' + args.chromosome + ':'), "Reading", args.cdnas_fasta)
        all_cds = read_fasta(args.cdnas_fasta, True)

        # align the variant coordinates to transcript, translate into the protein database
        print (('Chr ' + args.chromosome + ':'), 'Creating haplotype database.')
        haplo_results = process_haplotypes(all_transcripts, gene_haplo_df, all_cds, annotations_db, args.chromosome, args.haplo_id_prefix, args.force_rf, args.threads, args.min_freq, args.min_hap_count, args.ignore_UTR, args.skip_start_lost, (len(args.output_cdna_fasta) > 0))
        result_data = haplo_results[0]
        result_sequences = haplo_results[1]
        result_cdna = haplo_results[2]

        # store the result metadata        
        print (('Chr ' + args.chromosome + ':'), 'Storing the result metadata:', args.output_file)
        result_data.to_csv(args.output_file, sep='\t', header=True, index=False, compression='infer')
    
        # write the protein sequences into the fasta file
        output_fasta_file = gzip.open(args.output_fasta, 'wt') if args.output_fasta.endswith('.gz') else open(args.output_fasta, 'w')
        print (('Chr ' + args.chromosome + ':'), 'Writing FASTA file:', args.output_fasta)

        for i,seq in enumerate(result_sequences):
                accession = args.accession_prefix + '_' + hex(i)[2:]
                description = 'matching_proteins:' + ';'.join(seq['haplotypes']) + ' start:' + str(seq['start']) + ' reading_frame:' + ';'.join(seq['rfs'])

                output_fasta_file.write('>' + args.fasta_tag + '|' + accession + '|' + description + '\n')
                output_fasta_file.write(str(seq['sequence']) + '\n')

        output_fasta_file.close()

        # if requested, write the cDNA sequences into another fasta file
        if (len(args.output_cdna_fasta) > 0):
                output_cdna_fasta = gzip.open(args.output_cdna_fasta, 'wt') if args.output_cdna_fasta.endswith('.gz') else open(args.output_cdna_fasta, 'w')
                print (('Chr ' + args.chromosome + ':'), 'Writing cDNA FASTA file:', args.output_cdna_fasta)

                for i,seq in enumerate(result_cdna):
                       description = ';'.join(seq['haplotypes'])+ ' start:' + str(seq['start'])
                       output_cdna_fasta.write('>' + description + '\n')
                       output_cdna_fasta.write(str(seq['sequence']) + '\n')

                output_cdna_fasta.close()

        print (('Chr ' + args.chromosome + ':'), "Done.")
