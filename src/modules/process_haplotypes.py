import re
from multiprocessing import Pool
import pandas as pd
import bisect
from Bio.Seq import Seq
from modules.coordinates_toolbox import get_rna_position, get_rna_position_simple, check_start_change, get_affected_codons
from modules.common import KeyWrapper

result_columns = [  
    'TranscriptID',             # transcript stable ID from Ensembl
    'chromosome',               
    'transcript_biotype',       # transcript biotype annotation from Ensembl
    'HaplotypeID',              # artificial haplotype identifier -> matching in FASTA and result table
    'VCF_IDs',                  # IDs of corresponding lines in the VCF file
    'DNA_changes',              # list of changes in the DNA (could be in reversed order if gene is found on complementary strand)
    'allele_frequencies',       # list of minor allele frequencies for corresponding variants
    'cDNA_changes',             # list of changes mapped to the cDNA sequence
    'all_protein_changes',      # list of changes mapped to the protein sequence, including synonymous mutations
    'variant_types',            # type of variant (SAV, inframe-indel, synonymous, ...) for every change on the protein level
    'protein_changes',          # list of changes mapped to the protein sequence, excluding synonymous mutations
    'reading_frame',            # canonical reading frame for this transcript (-1 if unknown)
    'protein_prefix_length',    # length of the 5' UTR (in codons)    
    'start_missing',            # boolean - is the canonical annotation of a start codon present?
    'start_lost',               # boolean - does one of the changes cause loss of start codon?
    'splice_sites_affected',    # list of splice sites affected (1 means splice site between exon 1 annd 2, '-' if none)
#    'removed_DNA_changes',     # conflicting variants that have been removed - not used, transcripts where variants conflict are removed completely
    'occurrence_count',         # number of occurrences among the samples present in the VCF
    'frequency',                # frequency of occurrence among all samples present in the VCF
    'frequency_population',     # frequency of occurrence among the populations in the data set
    'frequency_superpopulation',# frequency of occurrence among the superpopulations in the data set
#    'samples'                   # samples containing this haplotype (in the format SAMPLE_ID:1 for maternal copy, SAMPLE_ID:2 for paternal copy) - not used currently, too large
]

# create dummy empty files in case of empty input
def empty_output(output_file, output_fasta):
    outfile = open(output_file, 'w')
    outfile.write('\t'.join(result_columns) + '\n')
    outfile.close()

    outfile = open(output_fasta, 'w')
    outfile.close()

# check if a protein variant maps to the UTR or not
# return boolean value
def check_protein_allele(change, start, stop):
    loc = int(change.split(':')[1].split('>')[1]) + start
    alt_len = len(change.split(':')[2])
    return (loc >= start) and (loc + alt_len <= stop)

def process_haplotypes(all_transcripts, genes_haplo_df, all_cdnas, annotations_db, chromosome, id_prefix, force_rf, threads, min_foo = -1, min_count = 0, ignore_UTR = True, skip_start_loss = True):
    result_data = []
    
    global process_transcript_haplotypes
    def process_transcript_haplotypes(transcript_feature):
        transcript_id = transcript_feature.id

        # Check if we have the cDNA sequence in the fasta
        if transcript_id not in all_cdnas:
            print('Transcript', transcript_id, 'not in cDNA database, skipping!')
            return []

        exons = [ exon for exon in annotations_db.children(transcript_feature, featuretype='exon', order_by='start') ]
        start_codons = [ sc for sc in annotations_db.children(transcript_feature, featuretype='start_codon', order_by='start') ]    # there should be only one, but just in case...
        stop_codons = [ sc for sc in annotations_db.children(transcript_feature, featuretype='stop_codon', order_by='start') ]      # not in use currently
        biotype = ""
        try:
            biotype = transcript_feature['transcript_biotype'][0]
        except:
            biotype = '-'

        # Some transcripts are classified as not coding -> start and stop codon positions are not given
        start_codon = None
        stop_codon = None

        if len(start_codons) > 0:
            start_codon = start_codons[0]
        if len(stop_codons) > 0:
            stop_codon = stop_codons[0]

        current_transcript = { 'ID': transcript_id, 'feature': transcript_feature, 'exons': exons, 'start_codon': start_codon, 'stop_codon': stop_codon, 'fasta_element': all_cdnas[transcript_id.split('.')[0]], 'biotype': biotype }
        transcript_haplotypes = genes_haplo_df[genes_haplo_df['TranscriptID'] == transcript_id]

        local_result_data = {}
        local_result_sequences = []      # way to avoid duplicate sequences -> access sequences by hash, aggregate haplotype IDs that correspond

        for index, row in transcript_haplotypes.iterrows():
            transcript_id = row['TranscriptID'].split('.')[0]

            # Check if any mutations are present
            if row['Changes'] == 'REF':
                continue

            cdna_sequence = current_transcript['fasta_element']['sequence'] # reference cDNA
            mutated_cdna = Seq(cdna_sequence)                               # cDNA to aggregate mutations
            sequence_length_diff = 0                                        # cummulative difference between the length of reference and haplotype -> to place mutations correctly in case of preceding indels

            # boolean - are we on a reverse strand?
            reverse_strand = current_transcript['feature'].strand == '-'

            cDNA_changes = []           # list of changes in the cDNA
            all_protein_changes = []    # list of changes in the protein sequence including synonymous mutations
            protein_changes = []        # list of changes in the protein sequence
            reading_frame = -1          # reading frame (0, 1 or 2), if known (inferred from the start codon position)
            reading_frame_ref = -1      # reading frame (0, 1 or 2) in the reference protein, if known (inferred from the start codon position) - can differ from the haplotype if there are indels in the UTR
            start_loc = 0               # location of the first nucleotide of the start codon with respect to the transcript start (0 if unknown)
            protein_start = 0           # length of the prefix (5'UTR) in protein
            protein_start_ref = 0       # length of the prefix (5'UTR) in the reference protein (can differ if there are indels in the UTR)
            spl_junctions_affected = [] # list of splicing junctions where a change takes place (identified by order, where 1 is the junction between the 1. and 2. exon), empty if none affected
            frameshifts = []            # boolean for every change whether it does or does not introduce a frameshift
            dna_var_types = []          # type of variant (SNP, indel, ...) for every change on the cDNA level
            prot_var_types = []         # type of variant (SAV, inframe-indel, ...) for every change on the protein level

            # Get the reading frame from the length between the start of the transcript and the start codon
            if (current_transcript['start_codon'] is not None):
                start_loc = get_rna_position_simple(transcript_id, current_transcript['start_codon'].start, current_transcript['exons'])
                if (reverse_strand):
                    start_loc = len(cdna_sequence) - start_loc - 3  

                reading_frame = start_loc % 3
                reading_frame_ref = reading_frame
                protein_start = int((start_loc - reading_frame) / 3)
                protein_start_ref = protein_start

            # Alternatively, use the stop codon in the same way, assume start at codon 0
            elif ((current_transcript['stop_codon'] is not None) and force_rf):
                stop_loc = get_rna_position_simple(transcript_id, current_transcript['stop_codon'].start, current_transcript['exons'])
                if (reverse_strand):
                    stop_loc = len(cdna_sequence) - stop_loc - 3  

                reading_frame = stop_loc % 3
                reading_frame_ref = reading_frame

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

                # get the general variant type
                if (mutation_intersects_intron is not None):
                    dna_var_types.append('splice')
                elif (ref_len == alt_len):
                    dna_var_types.append('SNP')
                else:
                    dna_var_types.append('indel')                

                # if we are on a reverse strand, we need to complement the reference and alternative sequence to match the cDNA
                # we also need to count the position from the end
                if reverse_strand:
                    ref_allele = ref_allele.reverse_complement()
                    alt_allele = alt_allele.reverse_complement()
                    rna_location = len(cdna_sequence) - rna_location - ref_len

                # check if the start codon gets shifted
                if (current_transcript['start_codon'] is not None) and (reading_frame > -1):
                    start_loc, reading_frame = check_start_change(start_loc, reading_frame, rna_location, ref_len, alt_len, force_rf)
                    if (start_loc == -1):
                        # if we wish to skip haplotypes where the start codon is lost, continue to another haplotype
                        if (skip_start_loss):
                            validity_check = False
                            break

                        dna_var_types[-1] = 'start_lost'    # replace the variant consequence
                        start_loc = 0
                        protein_start = 0
                        protein_start_ref = 0
                        reading_frame_ref = -1
                        start_lost = True
                    else:
                        protein_start = int((start_loc - reading_frame) / 3)

                # remember other derived attributes
                rna_locations.append(rna_location)
                ref_alleles.append(ref_allele)
                alt_alleles.append(alt_allele)
            
            if (not validity_check):
                continue

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
                rf_conseq = []      # protein consequence in each of the reading frames (e.g., SAV, frameshift, synonymous, etc.)
                is_synonymous = []
                
                ref_alleles_protein, protein_location_ref = get_affected_codons(cdna_sequence, rna_location, ref_len, reading_frame_ref, protein_start_ref)

                # need to readjust the position in the mutated sequence in case there were indels before
                rna_location += sequence_length_diff     

                alt_alleles_protein, protein_location_alt = get_affected_codons(mutated_cdna, rna_location, alt_len, reading_frame, protein_start)

                # store the change in protein as a string - only if there is a change (i.e. ignore synonymous variants) or a frameshift
                # loop through all possible reading frames
                for i,ref_allele_protein in enumerate(ref_alleles_protein):
                    alt_allele_protein = alt_alleles_protein[i]
                    loc_ref = protein_location_ref[i]
                    loc_alt = protein_location_alt[i]

                    if (dna_var_types[ch_idx] == 'splice'):
                        rf_conseq.append('splice_variant')
                    elif (dna_var_types[ch_idx] == 'start_lost'):
                        rf_conseq.append('start_lost')
                    elif (ref_allele_protein == alt_allele_protein):
                        rf_conseq.append('synonymous')
                    elif ('*' in ref_allele_protein) and not ('*' in alt_allele_protein):
                        rf_conseq.append('stop_lost') 
                    elif not ('*' in ref_allele_protein) and ('*' in alt_allele_protein):
                        rf_conseq.append('stop_gained')
                    elif (dna_var_types[ch_idx] == 'SNP'):
                        rf_conseq.append('SAV')
                    elif (dna_var_types[ch_idx] == 'indel'):
                        if (frameshifts[ch_idx]):
                            rf_conseq.append('frameshift')
                        else:
                            rf_conseq.append('inframe_indel')

                    is_synonymous.append(ref_allele_protein == alt_allele_protein)

                    protein_change = str(loc_ref) + ':' + ref_allele_protein + '>' + str(loc_alt) + ':' + alt_allele_protein
                    if (frameshifts[ch_idx]):
                        protein_change += "(+fs)"
                    elif (has_frameshift):
                        protein_change += "(fs)"
                        rf_conseq[-1] += '_after_fs'

                    rf_changes.append(protein_change)

                if not all(is_synonymous):
                    protein_changes.append("|".join(rf_changes))
                all_protein_changes.append("|".join(rf_changes))
                prot_var_types.append('|'.join(rf_conseq))

                has_frameshift = has_frameshift or frameshifts[ch_idx]
                sequence_length_diff += alt_len - ref_len
            
            spl_junctions_affected_str = ';'.join([str(x) for x in spl_junctions_affected])
            if (len(spl_junctions_affected_str) == 0):
                spl_junctions_affected_str = '-'

            haplotypeID = id_prefix + '_' + hex(index)[2:]

            # check the reading frame, if possible, and translate
            if (reading_frame > -1):
                protein_seq = mutated_cdna[reading_frame:].transcribe().translate()

                # since the reading frame is available, we know where the start codon is and it hasn't been lost by a mutation -> we can get rid of the UTRs
                # we don't filter UTR variants earlier, since changes in the start or stop codon redefine UTR regions in the transcript
                if (ignore_UTR):
                    # find the first stop codon after the start
                    first_stop = str(protein_seq).find('*', protein_start)

                    if (first_stop == -1):
                        first_stop = len(protein_seq)

                    # remove UTR variants
                    variant_filter = [ check_protein_allele(change, protein_start, first_stop) for change in all_protein_changes ]

                    all_vcf_IDs = [ vcf_id for idx,vcf_id in enumerate(all_vcf_IDs) if variant_filter[idx] ]
                    all_changes = [ ch for idx,ch in enumerate(all_changes) if variant_filter[idx] ]
                    all_AFs = [ af for idx,af in enumerate(all_AFs) if variant_filter[idx] ]
                    cDNA_changes = [ ch for idx,ch in enumerate(cDNA_changes) if variant_filter[idx] ]
                    all_protein_changes = [ ch for idx,ch in enumerate(all_protein_changes) if variant_filter[idx] ]
                    prot_var_types = [ t for idx,t in enumerate(prot_var_types) if variant_filter[idx] ]
                    protein_changes = [ ch for ch in protein_changes if check_protein_allele(ch, protein_start, first_stop) ]   # filter the only-non-synonymous variants separately as they are not indexed the same way

                    # skip this haplotype if no non-synonymous variants are left
                    if len (protein_changes) == 0:
                        continue

                    # cut away the 5' and 3' UTR sequences
                    protein_seq = protein_seq[protein_start:first_stop]
                    protein_start = 0

                # check if this haplotype is already in the results
                haplo_hash = str(hash(';'.join(all_vcf_IDs)))
                
                # if found, merge these two (increase the sample count and frequency)
                if haplo_hash in local_result_data:
                    local_result_data[haplo_hash][16] += row['Count']
                    local_result_data[haplo_hash][17] += row['Frequency']
                else:
                    local_result_data[haplo_hash] = [
                        row['TranscriptID'],
                        chromosome,
                        current_transcript['biotype'],
                        haplotypeID,
                        ';'.join(all_vcf_IDs),      # could be sorted in reverse order if on reverse strand
                        ";".join(all_changes),      
                        ";".join(all_AFs),
                        ';'.join(cDNA_changes),
                        ';'.join(all_protein_changes),    
                        ';'.join(prot_var_types),
                        ';'.join(protein_changes),
                        reading_frame,
                        protein_start,
                        current_transcript['start_codon'] is None,
                        start_lost,
                        spl_junctions_affected_str,
                        row['Count'],
                        row['Frequency'],
                        row['Frequency_population'],
                        row['Frequency_superpopulation'],
                    ]

                    # compute the sequence hash -> check if it already is in the list
                    # only store the sequence of the haplotype is not already stored
                    seq_hash = hash(str(protein_seq))
                    nearest_idx = bisect.bisect_left(KeyWrapper(local_result_sequences, key=lambda x: x['hash']), seq_hash)
                    if (len(local_result_sequences) > nearest_idx and local_result_sequences[nearest_idx]['hash'] == seq_hash):
                        local_result_sequences[nearest_idx]['haplotypes'].append(haplotypeID)
                        local_result_sequences[nearest_idx]['rfs'].append(str(reading_frame))
                    else:
                        local_result_sequences.insert(nearest_idx, {'hash': seq_hash, 'haplotypes': [haplotypeID], 'sequence': protein_seq, 'start': protein_start, 'rfs': [str(reading_frame)]})

            # unknown reading frame -> translate in all 3 reading frames
            # not possible to annotate UTRs -> keep everything
            elif (len(protein_changes) > 0):
                for rf in range(0,3):
                    protein_seq = mutated_cdna[rf:].transcribe().translate()

                    # compute the hash -> check if it already is in the list
                    seq_hash = hash(str(protein_seq))
                    nearest_idx = bisect.bisect_left(KeyWrapper(local_result_sequences, key=lambda x: x['hash']), seq_hash)
                    if (len(local_result_sequences) > nearest_idx and local_result_sequences[nearest_idx]['hash'] == seq_hash):
                        local_result_sequences[nearest_idx]['haplotypes'].append(haplotypeID)
                        local_result_sequences[nearest_idx]['rfs'].append(str(rf))
                    else:
                        local_result_sequences.insert(nearest_idx, {'hash': seq_hash, 'haplotypes': [haplotypeID], 'sequence': protein_seq, 'start': protein_start, 'rfs': [str(rf)]})

                # store result
                haplo_hash = hash(';'.join(all_vcf_IDs))

                local_result_data[haplo_hash] = [
                    row['TranscriptID'],
                    chromosome,
                    current_transcript['biotype'],
                    haplotypeID,
                    ';'.join(all_vcf_IDs),      # could be sorted in reverse order if on reverse strand
                    ";".join(all_changes),      
                    ";".join(all_AFs),
                    ';'.join(cDNA_changes),
                    ';'.join(all_protein_changes),  
                    ';'.join(prot_var_types),  
                    ';'.join(protein_changes),
                    reading_frame,
                    protein_start,
                    current_transcript['start_codon'] is None,
                    start_lost,
                    spl_junctions_affected_str,
                    row['Count'],
                    row['Frequency'],
                    row['Frequency_population'],
                    row['Frequency_superpopulation'],
                ]
                
        # filter haplotypes by frequency
        # filtering is done at the and since some haplotypes could have been merged
        result_table = []
        if (min_foo != -1):
            result_table = [ haplotype for haplotype in local_result_data.values() if (haplotype[-3] >= min_foo) ]
        else:
            result_table = [ haplotype for haplotype in local_result_data.values() if (haplotype[-4] >= min_count) ]

        included_haplotype_ids = [ haplotype[3] for haplotype in result_table ]

        # remove haplotypes below threshold from the FASTA header data 
        for seq in local_result_sequences:
            haplotypes_filter = [ (hap in included_haplotype_ids) for hap in seq['haplotypes'] ]    # list of booleans - are these haplotypes in the list after thresholding?
            seq['haplotypes'] = [ hap for hap_idx,hap in enumerate(seq['haplotypes']) if haplotypes_filter[hap_idx] ]
            seq['rfs'] = [ rf for hap_idx,rf in enumerate(seq['rfs']) if haplotypes_filter[hap_idx] ]

        # remove sequences where all matching haplotypes were filtered out
        local_result_sequences = [ seq for seq in local_result_sequences if len(seq['haplotypes']) > 0 ]

        return [result_table, local_result_sequences]

    #aggregated_results = list(map(process_transcript_haplotypes, all_transcripts))
    with Pool(threads) as p:
        aggregated_results = p.map(process_transcript_haplotypes, all_transcripts)

        result_data = []
        result_sequences = []

        for result_point in aggregated_results:
            if (len(result_point) == 2):
                result_data = result_data + list(result_point[0])
                result_sequences = result_sequences + result_point[1]

        result_df = pd.DataFrame(columns=result_columns, data=result_data)
        
        return [result_df, result_sequences]
