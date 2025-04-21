from datetime import datetime
import gzip
from numpy import ceil, floor
import pandas as pd
import bisect
from Bio.Seq import Seq
from modules.coordinates_toolbox import get_rna_position, get_rna_position_simple, check_start_change, get_affected_codons
from modules.common import KeyWrapper, check_vcf_df

result_columns = [  
    'transcriptID', 
    'chromosome',
    'transcript_biotype',
    'variantID',
    'vcfID',
    'DNA_change', 
    'cDNA_change', 
    'protein_change',
    'reading_frame', 
    'protein_prefix_length', 
    'start_missing',            # boolean - is the canonical annotation of a start codon present?
    'start_lost',               # boolean - does one of the changes cause loss of start codon?
    'splice_site_affected', 
]

# create dummy empty files in case of empty input
def empty_output(output_file, output_fasta):
    outfile = open(output_file, 'w')
    outfile.write('\t'.join(result_columns) + '\n')
    outfile.close()

    outfile = open(output_fasta, 'w')
    outfile.close()

# check if we gain a new stop codon
# return the location of the new start codon or -1
def check_start_gain(mutated_cdna, rna_location, alt_len):
    bpFrom = int(floor(rna_location / 3) * 3 )                  # first base of the first affected codon in RF 0                     
    bpTo = int(ceil((rna_location + alt_len - 2) / 3) * 3 + 2)  # last base of the last affected codon in RF 2

    # browse the potentially affected cDNA to search for an ATG codon
    for i in range(bpFrom, bpTo-3):
        codon = mutated_cdna[i:i+3]
        if str(codon) == 'ATG':
            return i

    return -1

def process_store_variants(all_transcripts, tmp_dir, log_file, all_cdnas, annotations_db, chromosome, fasta_tag, accession_prefix, force_rf, output_file, output_fasta, output_cdna_fasta):
    current_transcript = None
    result_data = []
    protein_sequence_list = []      # way to avoid duplicate sequences -> sort sequences by hash, aggregate variant IDs that correspond
    cdna_sequence_list = []      # way to avoid duplicate sequences -> sort sequences by hash, aggregate variant IDs that correspond

    for transcript in all_transcripts:
        transcript_id = transcript.id
        
        # load the according VCF file to Pandas
        vcf_df =  check_vcf_df(pd.read_csv(tmp_dir + '/' + transcript_id + '.tsv', sep='\t'))

        if (len(vcf_df) == 0):
            continue

        # Check if we have the cDNA sequence in the fasta
        if transcript_id not in all_cdnas:
            log_file.write('Transcript ' + transcript_id + ' does not have a reference cDNA sequence - skipping.\n')
            continue

        # store the annotation features of this transcript
        if (current_transcript is None or current_transcript['ID'] != transcript_id):
            transcript_feature = annotations_db[transcript_id]
            exons = [ exon for exon in annotations_db.children(transcript_feature, featuretype='exon', order_by='start') ]
            start_codons = [ sc for sc in annotations_db.children(transcript_feature, featuretype='start_codon', order_by='start') ]    # there should be only one, but just in case...
            stop_codons = [ sc for sc in annotations_db.children(transcript_feature, featuretype='stop_codon', order_by='start') ] 
            biotype = transcript_feature['transcript_biotype'][0]

            # start and stop codon positions are not given for some transcripts
            start_codon = None
            stop_codon = None

            if len(start_codons) > 0:
                start_codon = start_codons[0]
            if len(stop_codons) > 0:
                stop_codon = stop_codons[0]

            current_transcript = { 'ID': transcript_id, 'feature': transcript_feature, 'exons': exons, 'start_codon': start_codon, 'stop_codon': stop_codon, 'fasta_element': all_cdnas[transcript_id.split('.')[0]], 'biotype': biotype }
        
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

        # Alternatively, use the stop codon in the same way, assume start at codon 0
        elif ((current_transcript['stop_codon'] is not None) and force_rf):
            stop_loc = get_rna_position_simple(transcript_id, current_transcript['stop_codon'].start, current_transcript['exons'])
            if (reverse_strand):
                stop_loc = len(cdna_sequence) - stop_loc - 3  

            reading_frame = stop_loc % 3

        # iterate through rows of the VCF
        for index, vcf_row in vcf_df.iterrows():
            dna_location = int(vcf_row['POS'])

            ref_allele = ''
            if vcf_row['REF'] == '-':
                ref_allele = Seq('')
            else:
                ref_allele = Seq(vcf_row['REF'])

            alt_allele = ''
            if vcf_row['ALT'] == '-':
                alt_allele = Seq('')
            else:
                alt_allele = Seq(vcf_row['ALT'])

            var_ID = accession_prefix + '_' + transcript_id + '_' + hex(index)[2:]

            DNA_change = str(vcf_row['POS']) + ':' + vcf_row['REF'] + '>' + vcf_row['ALT']

            cDNA_change = ''                        # change in the cDNA
            protein_change = ''                     # change in the protein sequence
            spl_junction_affected = '-'             # splicing junction where a mutation takes place (identified by order, where 1 is the junction between the 1. and 2. exon), '-' if none affected
            mutated_cdna = Seq(cdna_sequence)       # cDNA sequence to mutate
            protein_start_variant = protein_start   # length of the UTR prefix in this protein variant (can differ if indel in 5' UTR)
            reading_frame_variant = reading_frame   # reading frame in this protein variant   
            start_lost = False                      # have we lost the canonical start codon?

            # compute the location in the RNA sequence
            # do any of the allele sequences intersect a splicing site? => truncate if so
            rna_location, ref_allele, ref_len, alt_allele, alt_len, spl_junction_affected = get_rna_position(transcript_id, dna_location, ref_allele, alt_allele, current_transcript['exons'])

            # if we are on a reverse strand, we need to complement the reference and alternative sequence to match the cDNA
            # we also need to count the position from the end
            if reverse_strand:
                ref_allele = ref_allele.reverse_complement()
                alt_allele = alt_allele.reverse_complement()
                rna_location = len(cdna_sequence) - rna_location - ref_len

            # check if what we expected to find is in fact in the cDNA
            if (str(ref_allele) != mutated_cdna[rna_location:rna_location+ref_len]):

                # if mismatched, check if shifting one base to the right will help
                dna_location += 1
                rna_location, ref_allele, ref_len, alt_allele, alt_len, spl_junction_affected = get_rna_position(transcript_id, dna_location, ref_allele, alt_allele, current_transcript['exons'])
                if reverse_strand:
                    rna_location = len(cdna_sequence) - rna_location - ref_len

                if (str(ref_allele) != mutated_cdna[rna_location:rna_location+ref_len]):

                    # if still mismatched, try one base to the left                    
                    dna_location -= 2
                    rna_location, ref_allele, ref_len, alt_allele, alt_len, spl_junction_affected = get_rna_position(transcript_id, dna_location, ref_allele, alt_allele, current_transcript['exons'])
                    if reverse_strand:
                        rna_location = len(cdna_sequence) - rna_location - ref_len

                    if (str(ref_allele) != mutated_cdna[rna_location:rna_location+ref_len]):
                        print('Ref allele not matching the cDNA sequence, skipping!')
                        log_file.write('[' + datetime.now().strftime('%X %x') + '] Ref allele not matching the cDNA sequence: ' + transcript_id + ' (reverse strand: ' + str(reverse_strand) + ') ID:' + vcf_row['ID'] + ' expected: ' + str(ref_allele) + ' found in cDNA: ' +  str(mutated_cdna[rna_location-10:rna_location]) + ' ' + str(mutated_cdna[rna_location:rna_location+ref_len]) + ' ' + str(mutated_cdna[rna_location+ref_len:rna_location+ref_len+10]) + '\n')
                        continue
            
            # apply the change to the cDNA
            mutated_cdna = mutated_cdna[:rna_location] + alt_allele + mutated_cdna[rna_location+ref_len:]

            cDNA_change = str(rna_location) + ':' + str(ref_allele) + '>' + str(alt_allele)

            # check if the start codon gets shifted
            if (current_transcript['start_codon'] is not None):
                new_start_loc, reading_frame_variant = check_start_change(start_loc, reading_frame, rna_location, ref_len, alt_len, force_rf)
                if (new_start_loc == -1):
                    protein_start_variant = 0
                    reading_frame_variant = -1
                    start_lost = True
                else:
                    protein_start_variant = int((new_start_loc - reading_frame) / 3)

            # remember reference allele in protein
            ref_alleles_protein = []        # reference residues directly affected (ignoring prossible frameshift), stored in a list for all three reading frames
            protein_location_ref = []       # location of these residues in the canonical protein (can be negative if in 5' UTR), creating a list as it can differ with reading frame

            ref_alleles_protein, protein_location_ref = get_affected_codons(cdna_sequence, rna_location, ref_len, reading_frame_variant, protein_start)
            
            alt_alleles_protein = []        # alternative allele in protein
            alt_alleles_protein, protein_location_alt = get_affected_codons(mutated_cdna, rna_location, alt_len, reading_frame_variant, protein_start_variant)

            # store the change in protein as a string - only if there is a change (i.e. ignore synonymous variants) or a frameshift
            # loop through all possible reading frames
            allele_changes = []
            for i,ref_allele_protein in enumerate(ref_alleles_protein):
                alt_allele_protein = alt_alleles_protein[i]
                loc_ref = protein_location_ref[i]
                loc_alt = protein_location_alt[i]

                change = str(loc_ref) + ':' + ref_allele_protein + '>' + str(loc_alt) + ':' + alt_allele_protein
                if (abs(ref_len - alt_len) % 3 > 0):
                    change += "(+fs)"

                allele_changes.append(change)

            # store the change in protein as a string - only if there is a change (i.e. ignore synonymous variants) or a frameshift
            protein_change = "|".join(allele_changes)

            # store result
            result_data.append([
                transcript_id,
                chromosome,
                current_transcript['biotype'],
                var_ID,
                vcf_row['ID'],
                DNA_change,
                cDNA_change,
                protein_change,
                reading_frame_variant,
                protein_start_variant,                
                current_transcript['start_codon'] is not None,
                start_lost,
                spl_junction_affected
            ])

            # check the reading frame, if possible, and translate
            if (reading_frame_variant > -1):
                protein_seq = mutated_cdna[reading_frame_variant:].transcribe().translate()

                # compute the hash -> check if it already is in the list
                seq_hash = hash(str(protein_seq))
                nearest_idx = bisect.bisect_left(KeyWrapper(protein_sequence_list, key=lambda x: x['hash']), seq_hash)
                if (len(protein_sequence_list) > nearest_idx and protein_sequence_list[nearest_idx]['hash'] == seq_hash):
                    protein_sequence_list[nearest_idx]['variants'].append(var_ID)
                    protein_sequence_list[nearest_idx]['rfs'].append(str(reading_frame_variant))
                else:
                    protein_sequence_list.insert(nearest_idx, {'hash': seq_hash, 'variants': [var_ID], 'sequence': protein_seq, 'start': protein_start_variant, 'rfs': [str(reading_frame_variant)]})


            # unknown reading frame -> translate in all 3 reading frames
            else:
                for rf in range(0,3):
                    protein_seq = mutated_cdna[rf:].transcribe().translate()

                    # compute the hash -> check if it already is in the list
                    seq_hash = hash(str(protein_seq))
                    nearest_idx = bisect.bisect_left(KeyWrapper(protein_sequence_list, key=lambda x: x['hash']), seq_hash)
                    if (len(protein_sequence_list) > nearest_idx and protein_sequence_list[nearest_idx]['hash'] == seq_hash):
                        protein_sequence_list[nearest_idx]['variants'].append(var_ID)
                        protein_sequence_list[nearest_idx]['rfs'].append(str(rf))
                    else:
                        protein_sequence_list.insert(nearest_idx, {'hash': seq_hash, 'variants': [var_ID], 'sequence': protein_seq, 'start': protein_start_variant, 'rfs': [str(rf)]})

            # If requested, store the cDNA sequence in the same way
            # compute the hash -> check if it already is in the list
            if (len(output_cdna_fasta) > 0):
                cdna_seq_hash = hash(str(mutated_cdna))
                nearest_idx = bisect.bisect_left(KeyWrapper(cdna_sequence_list, key=lambda x: x['hash']), cdna_seq_hash)
                if (len(cdna_sequence_list) > nearest_idx and cdna_sequence_list[nearest_idx]['hash'] == cdna_seq_hash):
                    cdna_sequence_list[nearest_idx]['variants'].append(var_ID)
                else:
                    cdna_sequence_list.insert(nearest_idx, {'hash': cdna_seq_hash, 'variants': [var_ID], 'sequence': mutated_cdna, 'start': max((reading_frame_variant + protein_start_variant * 3), -1)})

    # write the result table
    print ('Storing the result metadata:', output_file)
    result_df = pd.DataFrame(columns=result_columns, data=result_data)
    result_df.to_csv(output_file, sep='\t', header=True, index=False, compression='infer')

    # write the unique protein sequences into the fasta file
    output_fasta_file = gzip.open(output_fasta, 'wt') if (output_fasta.endswith('.gz')) else open(output_fasta, 'w')
    print ('Writing FASTA file:', output_fasta)

    for i,seq in enumerate(protein_sequence_list):
        accession = accession_prefix + '_chr' + chromosome + '_' + hex(i)[2:]
        description = 'matching_proteins:' + ';'.join(seq['variants']) + ' start:' + str(seq['start']) + ' reading_frame:' + ';'.join(seq['rfs'])

        output_fasta_file.write('>' + fasta_tag + '|' + accession + '|' + description + '\n')
        output_fasta_file.write(str(seq['sequence']) + '\n')

    output_fasta_file.close()
    
    if (len(output_cdna_fasta) > 0):
        outfile = gzip.open(output_cdna_fasta, 'wt') if (output_cdna_fasta.endswith('.gz')) else open(output_cdna_fasta, 'w')

        for seq in cdna_sequence_list:
            outfile.write('>' + ';'.join(seq['variants']) + ' start:' + str(seq['start']) + '\n')
            outfile.write(str(seq['sequence']) + '\n')

        outfile.close()