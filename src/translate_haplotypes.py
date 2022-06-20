from matplotlib.cbook import ls_mapper
from numpy import ceil, floor
import pandas as pd
import re
import argparse
import gffutils
from Bio.Seq import Seq
from common import read_fasta

# nucleotide_complement = { 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C' }

'''
Computes the position of the mutation in the RNA sequence. 
Checks whether the reference allele intersects a splice junction - truncates the sequences (both reference and alternative, if applicable) if so.
Special case: Mutation reaches over an intron into another exon. Probably covered but not tested.
'''
def get_rna_position(dna_location, ref_allele, alt_allele, exons):
    ref_len = len(ref_allele)
    alt_len = len(alt_allele)
    rna_location = 0
    found = False
    mutation_intersects_intron = None

    # find the corresponding exon - see how many nucleotides were there before
    for exon_idx,exon in enumerate(exons):
        # exon before the mutation -> remember the full length
        if (exon.end < dna_location):
            rna_location += (exon.end - exon.start + 1)
        
        # exon where the mutation happens -> remember position within
        # check for allele sequences ovrlapping borders of the exon
        elif (exon.start <= dna_location):
            rna_location += (dna_location - exon.start)
            found = True

            if (dna_location + ref_len > exon.end):
                remaining_length = exon.end - dna_location + 1
                mutation_intersects_intron = exon_idx + 1

                # check if the mutation does not reach into the next exon
                if exon_idx < (len(current_transcript['exons']) - 1) and (dna_location + ref_len > current_transcript['exons'][exon_idx+1].start):
                    next_exon = current_transcript['exons'][exon_idx+1]
                    start_again = next_exon.start - dna_location

                    ref_allele = ref_allele[:remaining_length] + ref_allele[start_again:]
                    ref_len = len(ref_allele)
                else:
                    ref_allele = ref_allele[:remaining_length]
                    ref_len = remaining_length

                # if there is an insertion that prolongs an exon, keep it, 
                # only truncate the alternative allele if the reference overlaps
                if (dna_location + alt_len > exon.end):
                    remaining_length = exon.end - dna_location + 1

                    # check if the mutation does not reach into the next exon
                    if exon_idx < (len(current_transcript['exons']) - 1) and (dna_location + alt_len > current_transcript['exons'][exon_idx+1].start):
                        next_exon = current_transcript['exons'][exon_idx+1]
                        start_again = next_exon.start - dna_location

                        alt_allele = alt_allele[:remaining_length] + alt_allele[start_again:]
                        alt_len = len(alt_allele)
                    else:
                        alt_allele = alt_allele[:remaining_length]
                        alt_len = remaining_length

            break

    if not found:
        raise Exception('DNA location ' + str(dna_location) + ' is not in an exon.')
    
    return rna_location, ref_allele, ref_len, alt_allele, alt_len, mutation_intersects_intron


def get_rna_position_simple(dna_location, exons):
    rna_location = 0
    found = False

    # find the corresponding exon - see how many nucleotides were there before
    for exon in exons:
        # exon before the location -> remember the full length
        if (exon.end < dna_location):
            rna_location += (exon.end - exon.start + 1)
        elif (exon.start <= dna_location):
            found = True
            rna_location += (dna_location - exon.start)
            break

    if not found:
        raise Exception('DNA location ' + str(dna_location) + ' is not in an exon.')

    return rna_location

parser = argparse.ArgumentParser(
    description="Processes a list of observed haplotypes in DNA -> computes location within transcripts and translated into protein haplotypes. Produces a FASTA file with all non-canonical protein sequences, and a CSV with all results and metadata.")

parser.add_argument("-db", dest="annotation_db", required=True,
                    help="DB file created by gffutils from GTF")

parser.add_argument("-i", dest="input_file", required=True,
                    help="input CSV file with gene haplotypes")

parser.add_argument("-cdna", dest="cdnas_fasta", required=True,
                    help="input cDNA fasta file")
                    
parser.add_argument("-tag", dest="fasta_tag", required=False,
                    help="tag for FASTa file entries", default='generic_enshap')

parser.add_argument("-acc_prefix", dest="accession_prefix", required=False,
                    help="tag for FASTa file entries", default='enshap_')

parser.add_argument("-output_csv", dest="output_file", required=True,
                    help="output CSV file")

parser.add_argument("-output_fasta", dest="output_fasta", required=True,
                    help="output FASTA file")

args = parser.parse_args()

'''
input_file = 'test/gene_haplo_chr22_v2.tsv'
output_file = 'test/haplo_chr22.tsv'
output_fasta = 'test/haplo_chr22.fa'
fasta_tag = 'generic_enshap'
annotations_file = 'test/chr22.db'
cdnas_fasta = 'test/total_cdnas.fa'
'''

output_fasta_file = open(args.output_fasta, 'w')

print ("Reading", args.input_file)
input_df = pd.read_csv(args.input_file, sep='\t')
print ("Reading", args.cdnas_fasta)
all_cdnas = read_fasta(args.cdnas_fasta)

# Load the annotations database
print ("Reading", args.annotations_db)
annotations_db = gffutils.FeatureDB(args.annotations_db)

print ('Processing haplotypes')
current_transcript = None
haplotypes_with_conflicting_mutations = 0

result_columns = [  
    'TranscriptID', 
    'HaplotypeID',
    'VCF_IDs', 
    'DNA_changes', 
    'allele_frequencies', 
    'cDNA_changes', 
    'protein_changes', 
    'reading_frame', 
    'protein_prefix_length', 
    'splice_sites_affected', 
    'removed_DNA_changes', 
    'occurrence_count', 
    'frequency', 
    'samples'
]

result_data = []

for index, row in input_df.iterrows():
    transcript_id = row['TranscriptID']

    # Check if any mutations are present
    if row['Changes'] == 'REF':
        continue

    # Check if we have the cDNA sequence in the fasta
    if transcript_id.split('.')[0] not in all_cdnas:
        continue

    # store the annotation features of this transcript
    if (current_transcript is None or current_transcript['ID'] != row['TranscriptID']):
        transcript_feature = annotations_db[transcript_id]
        exons = [ exon for exon in annotations_db.children(transcript_feature, featuretype='exon', order_by='start') ]
        start_codons = [ sc for sc in annotations_db.children(transcript_feature, featuretype='start_codon', order_by='start') ]    # there should be only one, but just in case...
        stop_codons = [ sc for sc in annotations_db.children(transcript_feature, featuretype='stop_codon', order_by='start') ]      # not in use currently

        # Some transcripts are classified as not coding -> start and stop codon positions are not given
        start_codon = None
        stop_codon = None

        if len(start_codons) > 0:
            start_codon = start_codons[0]
        if len(stop_codons) > 0:
            stop_codon = stop_codons[0]

        current_transcript = { 'ID': transcript_id, 'feature': transcript_feature, 'exons': exons, 'start_codon': start_codon, 'stop_codon': stop_codon, 'fasta_element': all_cdnas[transcript_id.split('.')[0]] }
    
    cdna_sequence = current_transcript['fasta_element']['sequence'] # reference cDNA
    mutated_cdna = cdna_sequence                                    # cDNA to aggregate mutations
    sequence_length_diff = 0                                        # cummulative difference between the length of reference and haplotype -> to place mutations correctly in case of preceding indels

    # boolean - are we on a reverse strand?
    reverse_strand = current_transcript['feature'].strand == '-'

    cDNA_changes = []           # list of changes in the cDNA
    protein_changes = []        # list of changes in the protein sequence
    reading_frame = -1          # reading frame (0, 1 or 2), if known (inferred from the start codon position)
    start_loc = 0               # location of the first nucleotide of the start codon with respect to the transcript start (0 if unknown)
    protein_start = 0           # length of the prefix in protein
    spl_junctions_affected = [] # list of splicing junctions where a mutation takes place (identified by order, where 1 is the junction between the 1. and 2. exon), empty if none affected
    #frameshifts = []            # boolean for every change whether it does or does not introduce a frameshift

    if (current_transcript['start_codon'] is not None):
        start_loc = get_rna_position_simple(current_transcript['start_codon'].start, current_transcript['exons'])
        reading_frame = start_loc % 3
        protein_start = int((start_loc - reading_frame) / 3)

    all_changes = row['Changes'].split(';')
    if (reverse_strand):
        all_changes.reverse()

    # iterate through mutations, construct mutated cdna
    for mutation in all_changes:
        ref_allele = Seq(re.split('\d+', mutation)[1].split('>')[0][1:])
        alt_allele = Seq(re.split('\d+', mutation)[1].split('>')[1])
        dna_location = int(re.split('[a-zA-Z\*]+', mutation)[0][:-1])

        # boolean - does this introduce a frameshift?
        # is_frameshift = (abs(len(ref_allele) - len(alt_allele)) % 3) != 0
        
        # compute the location in the RNA sequence
        # does any of the allele sequences intersect a splicing site? => truncate if so
        rna_location, ref_allele, ref_len, alt_allele, alt_len, mutation_intersects_intron = get_rna_position(dna_location, ref_allele, alt_allele, current_transcript['exons'])

        # if we are on a reverse strand, we need to complement the reference and alternative sequence to match the cDNA
        # we also need to count the position from the end
        if reverse_strand:
            ref_allele = ref_allele.reverse_complement()
            alt_allele = alt_allele.reverse_complement()
            rna_location = len(cdna_sequence) - rna_location - ref_len     
        
        # remember reference allele in protein
        ref_allele_protein = ""     # reference residues directly affected (ignoring prossible frameshift)
        protein_location = 0        # location of these residues in the canonical protein (can be negative if in 5' UTR)
        
        protein_location = int(floor((rna_location - max(reading_frame, 0)) / 3) - max(reading_frame, 0))           # if reading frame is unknown, assume 0
        bpFrom = int(floor((rna_location - max(reading_frame, 0)) / 3) * 3 + max(reading_frame, 0))
        bpTo = int(ceil((rna_location + len(ref_allele) - max(reading_frame, 0)) / 3) * 3 + max(reading_frame, 0))
        affected_codons = Seq(cdna_sequence[bpFrom:bpTo])
        ref_allele_protein = str(affected_codons.transcribe().translate())

        # store change in cDNA
        cDNA_changes.append(str(rna_location) + ':' + str(ref_allele) + '>' + str(alt_allele))

        #if not reverse_strand:
        # need to readjust the position in the mutated sequence in case there were indels before
        # no need on the reverse strand since we're changing from the end
        rna_location += sequence_length_diff     
           
        # check if what we expected to find is in fact in the cDNA
        if (str(ref_allele) != mutated_cdna[rna_location:rna_location+ref_len]):
            print('Ref allele not matching the cDNA sequence, skipping!')
            print(transcript_id + ' strand ' + current_transcript['feature'].strand + ' ' + mutation + ' cDNA: ' +  mutated_cdna[rna_location-10:rna_location] + ' ' + mutated_cdna[rna_location:rna_location+ref_len] + ' ' + mutated_cdna[rna_location+ref_len:rna_location+ref_len+10])
            cDNA_changes = cDNA_changes[:-1]    # remove the element that has already been added
            continue

        # is a splice junction affected? -> remember if so 
        if (mutation_intersects_intron is not None) and (mutation_intersects_intron not in spl_junctions_affected):
            spl_junctions_affected.append(mutation_intersects_intron)

        # apply the change to the cDNA
        mutated_cdna = mutated_cdna[:rna_location] + alt_allele + mutated_cdna[rna_location+ref_len:]
        sequence_length_diff += alt_len - ref_len

        # remember alternative allele in protein after the new location is computed
        # at this stage, it can only be done for the forward strand
        alt_allele_protein = ""
        
        bpFrom = int(floor((rna_location - max(reading_frame, 0)) / 3) * 3 + max(reading_frame, 0))
        bpTo = int(ceil((rna_location + len(alt_allele) - max(reading_frame, 0)) / 3) * 3 + max(reading_frame, 0))
        affected_codons = Seq(mutated_cdna[bpFrom:bpTo])
        alt_allele_protein = str(affected_codons.transcribe().translate())

        protein_change = str(protein_location) + ':' + ref_allele_protein + '>' + alt_allele_protein
        if ((sequence_length_diff % 3) > 0):
            protein_change += "(fs)"
            protein_changes.append(protein_change)
        elif (ref_allele_protein != alt_allele_protein):
            protein_changes.append(protein_change)


    cDNA_changes_str = ';'.join(cDNA_changes)

    protein_changes_str = ';'.join(protein_changes)
    #protein_changes_str = ""
    if (reading_frame == -1):
        protein_changes_str = 'UNKNOWN_REF'
    elif (len(protein_changes_str) == 0):
        protein_changes_str = 'REF'
    
    spl_junctions_affected_str = ';'.join(list(map(lambda x: str(x), spl_junctions_affected)))
    if (len(spl_junctions_affected_str) == 0):
        spl_junctions_affected_str = '-'

    haplotypeID = args.acc_prefix + hex(index)[2:]

    # store result
    result_data.append([
        row['TranscriptID'],
        haplotypeID,
        row['VCF_IDs'],
        row['Changes'],
        row['AlleleFrequencies'],
        cDNA_changes_str,
        protein_changes_str,
        reading_frame,
        protein_start,
        spl_junctions_affected_str,
        row['RemovedChanges'],
        row['Count'],
        row['Frequency'],
        row['Samples'],
    ])

    # check the reading frame, if possible, and translate
    if (reading_frame > -1 and protein_changes_str != 'REF'):
        protein_seq = Seq(mutated_cdna[reading_frame:]).transcribe().translate()

        output_fasta_file.write('>' + args.fasta_tag + '|' + haplotypeID + '|' + current_transcript['fasta_element']['description'] + ' start:' + str(protein_start) + '\n')
        output_fasta_file.write(str(protein_seq) + '\n')

    # unknown reading frame -> translate in all 3 reading frames
    elif (reading_frame == -1):
        for rf in range(0,3):
            protein_seq = Seq(mutated_cdna[rf:]).transcribe().translate()

            output_fasta_file.write('>' + args.fasta_tag + '|' + haplotypeID + '|' + current_transcript['fasta_element']['description'] + ' reading_frame:' + str(rf) + '\n')
            output_fasta_file.write(str(protein_seq) + '\n')


result_df = pd.DataFrame(columns=result_columns, data=result_data)
result_df.to_csv(args.output_file, sep='\t', header=True, index=False)

output_fasta_file.close()

print ("Done.")