import re
from numpy import ceil, floor
import pandas as pd
import bisect
from Bio.Seq import Seq
from coordinates_toolbox import get_rna_position, get_rna_position_simple
from common import KeyWrapper

result_columns = [  
    'TranscriptID', 
    'transcript_biotype',
    'HaplotypeID',
    'VCF_IDs', 
    'DNA_changes', 
    'allele_frequencies', 
    'cDNA_changes', 
    'all_protein_changes', 
    'protein_changes',
    'reading_frame', 
    'protein_prefix_length', 
    'splice_sites_affected', 
#    'removed_DNA_changes', 
    'occurrence_count', 
    'frequency', 
    'samples'
]

# create dummy empty files in case of empty input
def empty_output(output_file, output_fasta):
    outfile = open(output_file, 'w')
    outfile.write('\t'.join(result_columns) + '\n')
    outfile.close()

    outfile = open(output_fasta, 'w')
    outfile.close()

def process_store_haplotypes(genes_haplo_df, all_cdnas, annotations_db, chromosome, fasta_tag, id_prefix, accession_prefix, output_file, output_fasta):
    current_transcript = None
    result_data = []
    protein_sequence_list = []      # way to avoid duplicate sequences -> access sequences by hash, aggregate haplotype IDs that correspond
    
    for index, row in genes_haplo_df.iterrows():
        transcript_id = row['TranscriptID'].split('.')[0]

        # Check if any mutations are present
        if row['Changes'] == 'REF':
            continue

        # Check if we have the cDNA sequence in the fasta
        if transcript_id not in all_cdnas:
            continue

        # store the annotation features of this transcript
        if (current_transcript is None or current_transcript['ID'] != row['TranscriptID']):
            transcript_feature = annotations_db[transcript_id]
            exons = [ exon for exon in annotations_db.children(transcript_feature, featuretype='exon', order_by='start') ]
            start_codons = [ sc for sc in annotations_db.children(transcript_feature, featuretype='start_codon', order_by='start') ]    # there should be only one, but just in case...
            stop_codons = [ sc for sc in annotations_db.children(transcript_feature, featuretype='stop_codon', order_by='start') ]      # not in use currently
            biotype = transcript_feature['transcript_biotype'][0]

            # Some transcripts are classified as not coding -> start and stop codon positions are not given
            start_codon = None
            stop_codon = None

            if len(start_codons) > 0:
                start_codon = start_codons[0]
            if len(stop_codons) > 0:
                stop_codon = stop_codons[0]

            current_transcript = { 'ID': transcript_id, 'feature': transcript_feature, 'exons': exons, 'start_codon': start_codon, 'stop_codon': stop_codon, 'fasta_element': all_cdnas[transcript_id.split('.')[0]], 'biotype': biotype }
        
        cdna_sequence = current_transcript['fasta_element']['sequence'] # reference cDNA
        mutated_cdna = Seq(cdna_sequence)                               # cDNA to aggregate mutations
        sequence_length_diff = 0                                        # cummulative difference between the length of reference and haplotype -> to place mutations correctly in case of preceding indels

        # boolean - are we on a reverse strand?
        reverse_strand = current_transcript['feature'].strand == '-'

        cDNA_changes = []           # list of changes in the cDNA
        all_protein_changes = []    # list of changes in the protein sequence including synonymous mutations
        protein_changes = []        # list of changes in the protein sequence
        reading_frame = -1          # reading frame (0, 1 or 2), if known (inferred from the start codon position)
        start_loc = 0               # location of the first nucleotide of the start codon with respect to the transcript start (0 if unknown)
        protein_start = 0           # length of the prefix in protein
        spl_junctions_affected = [] # list of splicing junctions where a mutation takes place (identified by order, where 1 is the junction between the 1. and 2. exon), empty if none affected
        #frameshifts = []            # boolean for every change whether it does or does not introduce a frameshift

        # Get the reading frame from the length between the start of the transcript and the start codon
        if (current_transcript['start_codon'] is not None):
            start_loc = get_rna_position_simple(transcript_id, current_transcript['start_codon'].start, current_transcript['exons'])
            if (reverse_strand):
                start_loc = len(cdna_sequence) - start_loc - 3  

            reading_frame = start_loc % 3
            protein_start = int((start_loc - reading_frame) / 3)

        all_changes = row['Changes'].split(';')
        all_AFs = row['AlleleFrequencies'].split(';')
        all_vcf_IDs = row['VCF_IDs'].split(';')
        if (reverse_strand):                        # revert the order of mutations in the reverse strand, so that we add them from start to end and account for preceding indels
            all_AFs.reverse()
            all_changes.reverse()
            all_vcf_IDs.reverse()

        # iterate through mutations, construct mutated cdna
        for mutation in all_changes:
            ref_allele = Seq(re.split('\d+', mutation)[1].split('>')[0][1:])
            alt_allele = Seq(re.split('\d+', mutation)[1].split('>')[1])
            dna_location = int(re.split('[a-zA-Z\*]+', mutation)[0][:-1])

            # boolean - does this introduce a frameshift?
            # is_frameshift = (abs(len(ref_allele) - len(alt_allele)) % 3) != 0
            
            # compute the location in the RNA sequence
            # does any of the allele sequences intersect a splicing site? => truncate if so
            rna_location, ref_allele, ref_len, alt_allele, alt_len, mutation_intersects_intron = get_rna_position(transcript_id, dna_location, ref_allele, alt_allele, current_transcript['exons'])

            # if we are on a reverse strand, we need to complement the reference and alternative sequence to match the cDNA
            # we also need to count the position from the end
            if reverse_strand:
                ref_allele = ref_allele.reverse_complement()
                alt_allele = alt_allele.reverse_complement()
                rna_location = len(cdna_sequence) - rna_location - ref_len     
            
            # remember reference allele in protein
            ref_alleles_protein = []    # reference residues directly affected (ignoring prossible frameshift), stored in a list for all three reading frames
            protein_location_ref = []       # location of these residues in the canonical protein (can be negative if in 5' UTR), creating a list as it can differ with reading frame
            
            if (reading_frame == -1):
                for rf in range(3):                    
                    protein_location_ref.append(int(floor((rna_location - rf) / 3)))
            else:
                protein_location_ref = [int(floor((rna_location - reading_frame) / 3) -  protein_start)]

            bpFrom = int(floor((rna_location - max(reading_frame, 0)) / 3) * 3 + max(reading_frame, 0)) # if reading frame is unknown, assume 0 and add other reading frames later
            bpFrom = max(max(bpFrom, 0), reading_frame)                                                 # in case the beginning of the change is before the reading frame start

            bpTo = int(ceil((rna_location + len(ref_allele) - max(reading_frame, 0)) / 3) * 3 + max(reading_frame, 0))

            if (bpTo-bpFrom > 2): # make sure we have at least 1 codon covered
                affected_codons = Seq(cdna_sequence[bpFrom:bpTo])
                ref_alleles_protein = [str(affected_codons.transcribe().translate())]
            else:
                ref_alleles_protein = ['-']

            if reading_frame == -1:
                for rf in [1,2]:
                    bpFrom = int(floor((rna_location - rf) / 3) * 3 + rf) 
                    bpFrom = max(max(bpFrom, 0), rf)                                                    
                    bpTo = int(ceil((rna_location + len(ref_allele) - rf) / 3) * 3 + rf)

                    if (bpTo-bpFrom > 2): # make sure we have at least 1 codon covered
                        affected_codons = Seq(cdna_sequence[bpFrom:bpTo])
                        ref_alleles_protein.append(str(affected_codons.transcribe().translate()))
                    else:
                        ref_alleles_protein.append('-')


            # store change in cDNA
            cDNA_changes.append(str(rna_location) + ':' + str(ref_allele) + '>' + str(alt_allele))

            # need to readjust the position in the mutated sequence in case there were indels before
            rna_location += sequence_length_diff     
            
            # check if what we expected to find is in fact in the cDNA
            if (str(ref_allele) != mutated_cdna[rna_location:rna_location+ref_len]):
                print('Ref allele not matching the cDNA sequence, skipping!')
                print(transcript_id + ' strand ' + current_transcript['feature'].strand + ' ' + mutation + ' cDNA: ' +  mutated_cdna[rna_location-10:rna_location] + ' ' + mutated_cdna[rna_location:rna_location+ref_len] + ' ' + mutated_cdna[rna_location+ref_len:rna_location+ref_len+10])
                cDNA_changes = cDNA_changes[:-1]    # remove the element that has already been added
                continue

            # apply the change to the cDNA
            mutated_cdna = mutated_cdna[:rna_location] + alt_allele + mutated_cdna[rna_location+ref_len:]
            sequence_length_diff += alt_len - ref_len

            # remember alternative allele in protein after the new location is computed
            alt_alleles_protein = []
            protein_location_alt = []       # location of these residues in the canonical protein (can be negative if in 5' UTR), creating a list as it can differ with reading frame
            
            if (reading_frame == -1):
                for rf in range(3):                    
                    protein_location_alt.append(int(floor((rna_location - rf) / 3)))
            else:
                protein_location_alt = [int(floor((rna_location - reading_frame) / 3) -  protein_start)]
            
            bpFrom = int(floor((rna_location - max(reading_frame, 0)) / 3) * 3 + max(reading_frame, 0)) # if reading frame is unknown, assume 0 and add other reading frames later
            bpFrom = max(max(bpFrom, 0), reading_frame)                                                                     # in case the beginning of the change is before the reading frame start
            bpTo = int(ceil((rna_location + len(alt_allele) - max(reading_frame, 0)) / 3) * 3 + max(reading_frame, 0))

            if (bpTo-bpFrom > 2): # make sure we have at least 1 codon covered (i.e., the change doesn't fall before the reading frame start)
                affected_codons = mutated_cdna[bpFrom:bpTo]
                alt_alleles_protein = [str(affected_codons.transcribe().translate())]  
            else:
                alt_alleles_protein = ['-']

            if reading_frame == -1:
                for rf in [1,2]:
                    bpFrom = int(floor((rna_location - rf) / 3) * 3 + rf) 
                    bpFrom = max(max(bpFrom, 0), rf)                                                    
                    bpTo = int(ceil((rna_location + len(alt_allele) - rf) / 3) * 3 + rf)

                    if (bpTo-bpFrom > 2): # make sure we have at least 1 codon covered
                        affected_codons = mutated_cdna[bpFrom:bpTo]
                        alt_alleles_protein.append(str(affected_codons.transcribe().translate()))
                    else:
                        alt_alleles_protein.append('-')

            # store the change in protein as a string - only if there is a change (i.e. ignore synonymous variants) or a frameshift
            # loop through all possible reading frames
            allele_changes = []
            is_synonymous = []
            for i,ref_allele_protein in enumerate(ref_alleles_protein):
                alt_allele_protein = alt_alleles_protein[i]
                loc_ref = protein_location_ref[i]
                loc_alt = protein_location_alt[i]

                protein_change = str(loc_ref) + ':' + ref_allele_protein + '>' + str(loc_alt) + ':' + alt_allele_protein
                if (abs(len(ref_allele) - len(alt_allele)) % 3 > 0):
                    protein_change += "(+fs)"
                elif ((sequence_length_diff % 3) > 0):
                    protein_change += "(fs)"

                if (ref_allele_protein == alt_allele_protein):
                    is_synonymous.append(True)
                else:
                    is_synonymous.append(False)

                allele_changes.append(protein_change)

            if not all(is_synonymous):
                protein_changes.append("|".join(allele_changes))

            all_protein_changes.append("|".join(allele_changes))

            # is a splice junction affected? -> remember if so 
            if (mutation_intersects_intron is not None) and (mutation_intersects_intron not in spl_junctions_affected):
                spl_junctions_affected.append(mutation_intersects_intron)


        cDNA_changes_str = ';'.join(cDNA_changes)

        all_protein_changes_str = ';'.join(all_protein_changes)
        protein_changes_str = ';'.join(protein_changes)
        #protein_changes_str = ""
        if (len(protein_changes_str) == 0):
            protein_changes_str = 'REF'
        
        spl_junctions_affected_str = ';'.join(list(map(lambda x: str(x), spl_junctions_affected)))
        if (len(spl_junctions_affected_str) == 0):
            spl_junctions_affected_str = '-'

        haplotypeID = id_prefix + '_' + hex(index)[2:]

        # store result
        result_data.append([
            row['TranscriptID'],
            current_transcript['biotype'],
            haplotypeID,
            ';'.join(all_vcf_IDs),      # could be sorted in reverse order if on reverse strand
            ";".join(all_changes),      
            ";".join(all_AFs),
            cDNA_changes_str,
            all_protein_changes_str,
            protein_changes_str,
            reading_frame,
            protein_start,
            spl_junctions_affected_str,
#            row['RemovedChanges'],
            row['Count'],
            row['Frequency'],
            row['Samples'],
        ])

        # check the reading frame, if possible, and translate
        if (reading_frame > -1 and protein_changes_str != 'REF'):
            protein_seq = mutated_cdna[reading_frame:].transcribe().translate()

            # compute the hash -> check if it already is in the list
            seq_hash = hash(str(protein_seq))
            nearest_idx = bisect.bisect_left(KeyWrapper(protein_sequence_list, key=lambda x: x['hash']), seq_hash)
            if (len(protein_sequence_list) > nearest_idx and protein_sequence_list[nearest_idx]['hash'] == seq_hash):
                protein_sequence_list[nearest_idx]['haplotypes'].append(haplotypeID)
                protein_sequence_list[nearest_idx]['rfs'].append(str(reading_frame))
            else:
                protein_sequence_list.insert(nearest_idx, {'hash': seq_hash, 'haplotypes': [haplotypeID], 'sequence': protein_seq, 'start': protein_start, 'rfs': [str(reading_frame)]})

            #output_fasta_file.write('>' + fasta_tag + '|' + haplotypeID + '|' + current_transcript['fasta_element']['description'] + ' start:' + str(protein_start) + '\n')
            #output_fasta_file.write(str(protein_seq) + '\n')

        # unknown reading frame -> translate in all 3 reading frames
        elif (protein_changes_str != 'REF'):
            for rf in range(0,3):
                protein_seq = mutated_cdna[rf:].transcribe().translate()

                # compute the hash -> check if it already is in the list
                seq_hash = hash(str(protein_seq))
                nearest_idx = bisect.bisect_left(KeyWrapper(protein_sequence_list, key=lambda x: x['hash']), seq_hash)
                if (len(protein_sequence_list) > nearest_idx and protein_sequence_list[nearest_idx]['hash'] == seq_hash):
                    protein_sequence_list[nearest_idx]['haplotypes'].append(haplotypeID)
                    protein_sequence_list[nearest_idx]['rfs'].append(str(rf))
                else:
                    protein_sequence_list.insert(nearest_idx, {'hash': seq_hash, 'haplotypes': [haplotypeID], 'sequence': protein_seq, 'start': protein_start, 'rfs': [str(rf)]})

                #output_fasta_file.write('>' + fasta_tag + '|' + haplotypeID + '|' + current_transcript['fasta_element']['description'] + ' reading_frame:' + str(rf) + '\n')
                #output_fasta_file.write(str(protein_seq) + '\n')


    print ('Storing the result metadata:', output_file)
    result_df = pd.DataFrame(columns=result_columns, data=result_data)
    result_df.to_csv(output_file, sep='\t', header=True, index=False)
    
    # write the unique protein sequences into the fasta file
    output_fasta_file = open(output_fasta, 'w')
    print ('Writing FASTA file:', output_fasta)

    for i,seq in enumerate(protein_sequence_list):
        accession = accession_prefix + '_' + hex(i)[2:]
        description = 'matching_proteins:' + ';'.join(seq['haplotypes']) + ' start:' + str(seq['start']) + ' reading_frame:' + ';'.join(seq['rfs'])

        output_fasta_file.write('>' + fasta_tag + '|' + accession + '|' + description + '\n')
        output_fasta_file.write(str(seq['sequence']) + '\n')

    output_fasta_file.close()

    print ("Done.")
