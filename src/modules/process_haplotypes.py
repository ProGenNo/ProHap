import re
from numpy import ceil, floor
import pandas as pd
import bisect
from Bio.Seq import Seq
from modules.coordinates_toolbox import get_rna_position, get_rna_position_simple, check_start_change, get_affected_codons
from modules.common import KeyWrapper

result_columns = [  
    'TranscriptID',             # transcript stable ID from Ensembl
    'transcript_biotype',       # transcript biotype annotation from Ensembl
    'HaplotypeID',              # artificial haplotype identifier -> matching in FASTA and result table
    'VCF_IDs',                  # IDs of corresponding lines in the VCF file
    'DNA_changes',              # list of changes in the DNA (could be in reversed order if gene is found on complementary strand)
    'allele_frequencies',       # list of minor allele frequencies for corresponding variants
    'cDNA_changes',             # list of changes mapped to the cDNA sequence
    'all_protein_changes',      # list of changes mapped to the protein sequence, including synonymous mutations
    'protein_changes',          # list of changes mapped to the protein sequence, excluding synonymous mutations
    'reading_frame',            # canonical reading frame for this transcript (-1 if unknown)
    'protein_prefix_length',    # length of the 5' UTR (in codons)    
    'start_lost',               # boolean - does one of the changes cause loss of start codon?
    'splice_sites_affected',    # list of splice sites affected (1 means splice site between exon 1 annd 2, '-' if none)
#    'removed_DNA_changes', 
    'occurrence_count',         # number of occurrences among the samples present in the VCF
    'frequency',                # frequency of occurrence among the samples present in the VCF
    'samples'                   # samples containing this haplotype (in the format SAMPLE_ID:1 for maternal copy, SAMPLE_ID:2 for paternal copy)
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
        spl_junctions_affected = [] # list of splicing junctions where a change takes place (identified by order, where 1 is the junction between the 1. and 2. exon), empty if none affected
        frameshifts = []            # boolean for every change whether it does or does not introduce a frameshift

        # Get the reading frame from the length between the start of the transcript and the start codon
        if (current_transcript['start_codon'] is not None):
            start_loc = get_rna_position_simple(transcript_id, current_transcript['start_codon'].start, current_transcript['exons'])
            if (reverse_strand):
                start_loc = len(cdna_sequence) - start_loc - 3  

            reading_frame = start_loc % 3
            protein_start = int((start_loc - reading_frame) / 3)

        # Alternatively, use the stop codon in the same way, assume start at codon 0
        if (current_transcript['stop_codon'] is not None):
            stop_loc = get_rna_position_simple(transcript_id, current_transcript['stop_codon'].start, current_transcript['exons'])
            if (reverse_strand):
                stop_loc = len(cdna_sequence) - stop_loc - 3  

            reading_frame = stop_loc % 3

        all_changes = row['Changes'].split(';')
        all_AFs = row['AlleleFrequencies'].split(';')
        all_vcf_IDs = row['VCF_IDs'].split(';')
        if (reverse_strand):                        # revert the order of mutations in the reverse strand, so that we add them from start to end and account for preceding indels
            all_AFs.reverse()
            all_changes.reverse()
            all_vcf_IDs.reverse()

        ref_alleles = []
        alt_alleles = []
        rna_locations = []
        start_lost = False
        validity_check = True

        # iterate through changes first time, check if start codon is lost, 
        # adjust the position of the start codon if shifted (i.e. inframe indel in 5' UTR)
        # remember alleles, locations on DNA and RNA
        for change in all_changes:
            ref_allele = re.split('\d+', change)[1].split('>')[0][1:]
            alt_allele = re.split('\d+', change)[1].split('>')[1]

            # in case the allele is fully deleted
            if ref_allele == '-':
                ref_allele = Seq('')
            else:
                ref_allele = Seq(ref_allele)

            if alt_allele == '-':
                alt_allele = Seq('')
            else:
                alt_allele = Seq(alt_allele)

            dna_location = int(re.split('[a-zA-Z\*]+', change)[0][:-1])
            
            # compute the location in the RNA sequence
            # does any of the allele sequences intersect a splicing site? => truncate if so
            rna_location, ref_allele, ref_len, alt_allele, alt_len, mutation_intersects_intron = get_rna_position(transcript_id, dna_location, ref_allele, alt_allele, current_transcript['exons'])

            # is a splice junction affected? -> remember if so 
            if (mutation_intersects_intron is not None) and (mutation_intersects_intron not in spl_junctions_affected):
                spl_junctions_affected.append(int(mutation_intersects_intron))

            # boolean - does this introduce a frameshift?
            frameshifts.append((abs(ref_len - alt_len) % 3) != 0)

            # if we are on a reverse strand, we need to complement the reference and alternative sequence to match the cDNA
            # we also need to count the position from the end
            if reverse_strand:
                ref_allele = ref_allele.reverse_complement()
                alt_allele = alt_allele.reverse_complement()
                rna_location = len(cdna_sequence) - rna_location - ref_len

            # check if the start codon gets shifted
            if (current_transcript['start_codon'] is not None) and (reading_frame > -1):
                start_loc = check_start_change(start_loc, rna_location, ref_len, alt_len)
                if (start_loc == -1):
                    start_loc = 0
                    protein_start = 0
                    reading_frame = -1
                    start_lost = True
                else:
                    protein_start = int((start_loc - reading_frame) / 3)

            # remember other derived attributes
            rna_locations.append(rna_location)
            ref_alleles.append(ref_allele)
            alt_alleles.append(alt_allele)

        # iterate through mutations second time, construct mutated cdna
        for ch_idx,change in enumerate(all_changes):
            ref_allele = ref_alleles[ch_idx]
            alt_allele = alt_alleles[ch_idx]
            ref_len = len(ref_allele)
            alt_len = len(alt_allele)

            rna_location = rna_locations[ch_idx]

            # store change in cDNA
            cDNA_changes.append(str(rna_location) + ':' + str(ref_allele) + '>' + str(alt_allele))

            # need to readjust the position in the mutated sequence in case there were indels before
            rna_location += sequence_length_diff     
            
            # check if what we expected to find is in fact in the cDNA
            if (str(ref_allele) != mutated_cdna[rna_location:rna_location+ref_len]):
                print('Ref allele not matching the cDNA sequence, skipping haplotype!')
                print(transcript_id + ' strand ' + current_transcript['feature'].strand + ' ' + change + ' cDNA: ' +  mutated_cdna[rna_location-10:rna_location] + ' ' + mutated_cdna[rna_location:rna_location+ref_len] + ' ' + mutated_cdna[rna_location+ref_len:rna_location+ref_len+10])
                cDNA_changes = cDNA_changes[:-1]    # remove the element that has already been added
                validity_check = False
                break

            # apply the change to the cDNA
            mutated_cdna = mutated_cdna[:rna_location] + alt_allele + mutated_cdna[rna_location+ref_len:]
            sequence_length_diff += alt_len - ref_len

        # skip in case of misaligned change
        if (not validity_check):
            continue

        # translate all the individual changes only after all the changes have been added to the mutated cDNA -> cover cases where multiple mutations affect the same codon(s)
        cdna_sequence = Seq(cdna_sequence)
        has_frameshift = False
        sequence_length_diff = 0

        # iterate third time, remember individual changes on the protein level
        for ch_idx,change in enumerate(all_changes):            
            ref_allele = ref_alleles[ch_idx]
            alt_allele = alt_alleles[ch_idx]
            ref_len = len(ref_allele)
            alt_len = len(alt_allele)

            rna_location = rna_locations[ch_idx]

            rf_changes = []     # residues directly affected (ignoring prossible frameshift), stored in a list for all three reading frames
            is_synonymous = []
            
            ref_alleles_protein, protein_location_ref = get_affected_codons(cdna_sequence, rna_location, ref_len, reading_frame, protein_start)

            # need to readjust the position in the mutated sequence in case there were indels before
            rna_location += sequence_length_diff     

            alt_alleles_protein, protein_location_alt = get_affected_codons(mutated_cdna, rna_location, alt_len, reading_frame, protein_start)

            # store the change in protein as a string - only if there is a change (i.e. ignore synonymous variants) or a frameshift
            # loop through all possible reading frames
            for i,ref_allele_protein in enumerate(ref_alleles_protein):
                alt_allele_protein = alt_alleles_protein[i]
                loc_ref = protein_location_ref[i]
                loc_alt = protein_location_alt[i]

                is_synonymous.append(ref_allele_protein == alt_allele_protein)

                protein_change = str(loc_ref) + ':' + ref_allele_protein + '>' + str(loc_alt) + ':' + alt_allele_protein
                if (frameshifts[i]):
                    protein_change += "(+fs)"
                elif (has_frameshift):
                    protein_change += "(fs)"

                rf_changes.append(protein_change)

            if not all(is_synonymous):
                protein_changes.append("|".join(rf_changes))
            all_protein_changes.append("|".join(rf_changes))

            has_frameshift = has_frameshift or frameshifts[i]
            sequence_length_diff += alt_len - ref_len

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
            start_lost,
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
