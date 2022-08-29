import re
from numpy import ceil, floor
import pandas as pd
import bisect
from Bio.Seq import Seq
from coordinates_toolbox import get_rna_position, get_rna_position_simple
from common import KeyWrapper

result_columns = [  
    'transcriptID', 
    'chromosome',
    'transcript_biotype',
    'variantID',
    'DNA_change', 
    'allele_frequency', 
    'cDNA_change', 
    'protein_change',
    'reading_frame', 
    'protein_prefix_length', 
    'splice_site_affected', 
]

def process_store_variants(all_transcripts, tmp_dir, log_file, all_cdnas, annotations_db, chromosome, fasta_tag, accession_prefix, output_file, output_fasta):
    current_transcript = None
    result_data = []
    protein_sequence_list = []      # way to avoid duplicate sequences -> access sequences by hash, aggregate variant IDs that correspond

    for transcript in all_transcripts:
        transcript_id = transcript.id

        # Check if we have the cDNA sequence in the fasta
        if transcript_id not in all_cdnas:
            log_file.write('Transcript ' + transcript_id + ' does not have a reference cDNA sequence - skipping.\n')
            continue

        # store the annotation features of this transcript
        if (current_transcript is None or current_transcript['ID'] != transcript_id):
            transcript_feature = annotations_db[transcript_id]
            exons = [ exon for exon in annotations_db.children(transcript_feature, featuretype='exon', order_by='start') ]
            start_codons = [ sc for sc in annotations_db.children(transcript_feature, featuretype='start_codon', order_by='start') ]    # there should be only one, but just in case...
            stop_codons = [ sc for sc in annotations_db.children(transcript_feature, featuretype='stop_codon', order_by='start') ]      # not in use currently
            biotype = transcript_feature['transcript_biotype'][0]

            # start and stop codon positions are not given for some transcripts
            start_codon = None
            stop_codon = None

            if len(start_codons) > 0:
                start_codon = start_codons[0]
            if len(stop_codons) > 0:
                stop_codon = stop_codons[0]

            current_transcript = { 'ID': transcript_id, 'feature': transcript_feature, 'exons': exons, 'start_codon': start_codon, 'stop_codon': stop_codon, 'fasta_element': all_cdnas[transcript_id.split('.')[0]], 'biotype': biotype }
        
        # load the according VCF file to Pandas
        vcf_df = pd.read_csv(tmp_dir + '/' + transcript_id + '.tsv', sep='\t')

        cdna_sequence = current_transcript['fasta_element']['sequence']  # reference cDNA        
        reverse_strand = current_transcript['feature'].strand == '-'     # boolean - are we on a reverse strand?

        reading_frame = -1          # reading frame (0, 1 or 2), if known (inferred from the start codon position), -1 if unknown
        start_loc = 0               # location of the first nucleotide of the start codon with respect to the transcript start (0 if unknown)
        protein_start = 0           # length of the prefix in protein     

        # Get the reading frame from the length between the start of the transcript and the start codon
        if (current_transcript['start_codon'] is not None):
            start_loc = get_rna_position_simple(transcript_id, current_transcript['start_codon'].start, current_transcript['exons'])
            if (reverse_strand):
                start_loc = len(cdna_sequence) - start_loc - 3  

            reading_frame = start_loc % 3
            protein_start = int((start_loc - reading_frame) / 3)

        # iterate through rows of the VCF
        for index, vcf_row in vcf_df.iterrows():
            dna_location = int(vcf_row['POS'])
            ref_allele = Seq(vcf_row['REF'])
            alt_allele = Seq(vcf_row['ALT'])

            DNA_change = vcf_row['POS'] + ':' + vcf_row['REF'] + '>' + vcf_row['ALT']

            cDNA_change = ''                    # change in the cDNA
            protein_change = ''                 # change in the protein sequence
            spl_junction_affected = '-'         # splicing junction where a mutation takes place (identified by order, where 1 is the junction between the 1. and 2. exon), '-' if none affected
            mutated_cdna = Seq(cdna_sequence)   # cDNA sequence to mutate

            # compute the location in the RNA sequence
            # does any of the allele sequences intersect a splicing site? => truncate if so
            rna_location, ref_allele, ref_len, alt_allele, alt_len, spl_junction_affected = get_rna_position(transcript_id, dna_location, ref_allele, alt_allele, current_transcript['exons'])

            # if we are on a reverse strand, we need to complement the reference and alternative sequence to match the cDNA
            # we also need to count the position from the end
            if reverse_strand:
                ref_allele = ref_allele.reverse_complement()
                alt_allele = alt_allele.reverse_complement()
                rna_location = len(cdna_sequence) - rna_location - ref_len

            # remember reference allele in protein
            ref_allele_protein = ""     # reference residues directly affected (ignoring prossible frameshift)
            protein_location = 0        # location of these residues in the canonical protein (can be negative if in 5' UTR)
            
            if (reading_frame == -1):
                protein_location = int(floor((rna_location - max(reading_frame, 0)) / 3))               # if reading frame is unknown, assume 0
            else:
                protein_location = int(floor((rna_location - max(reading_frame, 0)) / 3) -  protein_start)

            bpFrom = int(floor((rna_location - max(reading_frame, 0)) / 3) * 3 + max(reading_frame, 0)) # if reading frame is unknown, assume 0
            bpFrom = max(bpFrom, 0)                                                                     # in case the beginning of the change is before the reading frame start

            bpTo = int(ceil((rna_location + len(ref_allele) - max(reading_frame, 0)) / 3) * 3 + max(reading_frame, 0))

            if (bpTo >= 2): # make sure we have at least 1 codon covered
                affected_codons = Seq(cdna_sequence[bpFrom:bpTo])
                ref_allele_protein = str(affected_codons.transcribe().translate())

            cDNA_change = str(rna_location) + ':' + str(ref_allele) + '>' + str(alt_allele)

            # check if what we expected to find is in fact in the cDNA
            if (str(ref_allele) != mutated_cdna[rna_location:rna_location+ref_len]):
                print('Ref allele not matching the cDNA sequence, skipping!')
                log_file.write('Ref allele not matching the cDNA sequence: ' + transcript_id + ' strand ' + current_transcript['feature'].strand + ' ' + DNA_change + ' cDNA: ' +  mutated_cdna[rna_location-10:rna_location] + ' ' + mutated_cdna[rna_location:rna_location+ref_len] + ' ' + mutated_cdna[rna_location+ref_len:rna_location+ref_len+10] + '\n')
                continue

            # apply the change to the cDNA
            mutated_cdna = mutated_cdna[:rna_location] + alt_allele + mutated_cdna[rna_location+ref_len:]
            sequence_length_diff += alt_len - ref_len

            # remember alternative allele in protein after the new location is computed
            # at this stage, it can only be done for the forward strand
            alt_allele_protein = ""
            
            bpFrom = int(floor((rna_location - max(reading_frame, 0)) / 3) * 3 + max(reading_frame, 0)) # if reading frame is unknown, assume 0
            bpFrom = max(bpFrom, 0)                                                                     # in case the beginning of the change is before the reading frame start
            bpTo = int(ceil((rna_location + len(alt_allele) - max(reading_frame, 0)) / 3) * 3 + max(reading_frame, 0))

            if (bpTo >= 2): # make sure we have at least 1 codon covered (i.e., the change doesn't fall before the reading frame start)
                affected_codons = mutated_cdna[bpFrom:bpTo]
                alt_allele_protein = str(affected_codons.transcribe().translate())

            # store the change in protein as a string - only if there is a change (i.e. ignore synonymous variants) or a frameshift
            protein_change = str(protein_location) + ':' + ref_allele_protein + '>' + alt_allele_protein
            if (abs(len(ref_allele) - len(alt_allele)) % 3 > 0):
                protein_change += "(+fs)"

            # store result
            result_data.append([
                transcript_id,
                chromosome,
                current_transcript['biotype'],
                vcf_row['ID'],
                DNA_change,
                cDNA_change,
                protein_change,
                reading_frame,
                protein_start,
                spl_junction_affected
            ])

            # check the reading frame, if possible, and translate
            if (reading_frame > -1):
                protein_seq = mutated_cdna[reading_frame:].transcribe().translate()

                # compute the hash -> check if it already is in the list
                seq_hash = hash(str(protein_seq))
                nearest_idx = bisect.bisect_left(KeyWrapper(protein_sequence_list, key=lambda x: x['hash']), seq_hash)
                if (len(protein_sequence_list) > nearest_idx and protein_sequence_list[nearest_idx]['hash'] == seq_hash):
                    protein_sequence_list[nearest_idx]['variants'].append(vcf_row['ID'])
                    protein_sequence_list[nearest_idx]['rfs'].append(str(reading_frame))
                else:
                    protein_sequence_list.insert(nearest_idx, {'hash': seq_hash, 'variants': [vcf_row['ID']], 'sequence': protein_seq, 'start': protein_start, 'rfs': [str(reading_frame)]})


            # unknown reading frame -> translate in all 3 reading frames
            else:
                for rf in range(0,3):
                    protein_seq = mutated_cdna[rf:].transcribe().translate()

                    # compute the hash -> check if it already is in the list
                    seq_hash = hash(str(protein_seq))
                    nearest_idx = bisect.bisect_left(KeyWrapper(protein_sequence_list, key=lambda x: x['hash']), seq_hash)
                    if (len(protein_sequence_list) > nearest_idx and protein_sequence_list[nearest_idx]['hash'] == seq_hash):
                        protein_sequence_list[nearest_idx]['variants'].append(vcf_row['ID'])
                        protein_sequence_list[nearest_idx]['rfs'].append(str(rf))
                    else:
                        protein_sequence_list.insert(nearest_idx, {'hash': seq_hash, 'variants': [vcf_row['ID']], 'sequence': protein_seq, 'start': protein_start, 'rfs': [str(rf)]})

    # write the result table
    print ('Storing the result metadata:', output_file)
    result_df = pd.DataFrame(columns=result_columns, data=result_data)
    result_df.to_csv(output_file, sep='\t', header=True, index=False)

    # write the unique protein sequences into the fasta file
    output_fasta_file = open(output_fasta, 'w')
    print ('Writing FASTA file:', output_fasta)

    for i,seq in enumerate(protein_sequence_list):
        accession = accession_prefix + '_' + hex(i)[2:]
        description = 'matching_proteins:' + ';'.join(seq['variants']) + ' start:' + str(seq['start']) + ' reading_frame:' + ';'.join(seq['rfs'])

        output_fasta_file.write('>' + fasta_tag + '|' + accession + '|' + description + '\n')
        output_fasta_file.write(str(seq['sequence']) + '\n')

    output_fasta_file.close()




            

            


