import argparse
import pandas as pd
import re
import bisect
import gffutils
from datetime import datetime
from tqdm import tqdm
from multiprocessing import Pool
from common import read_fasta

parser = argparse.ArgumentParser(
	description='Reads the peptide report file, creates a new file with information about covered SNPs for each peptide and protein')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input CSV file")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

parser.add_argument("-hap_tsv", dest="haplo_db", required=False,
                    help="haplotypes tab-separated file (optional)", default=None)

parser.add_argument("-var_tsv", dest="var_db", required=False,
                    help="variants tab-separated file (optional)", default=None)

parser.add_argument("-hap_prefix", dest="haplo_prefix", required=False,
                    help="prefix for haplotype protein ID (default: 'haplo_')", default='haplo_')

parser.add_argument("-var_prefix", dest="var_prefix", required=False,
                    help="prefix for variant protein ID (default: 'var_')", default='var_')

parser.add_argument("-manual_prefix", dest="manual_prefix", required=False,
                    help="prefix for manually added protein sequences (default: 'man_')", default='man_')

parser.add_argument("-ens_annot", dest="annotation_db", required=True,
                    help="DB file created by gffutils from Ensembl GTF")

parser.add_argument("-t", dest="threads", type=int, required=True,
                    help="maximum number of threads")

parser.add_argument("-ref_fa", dest="ref_fasta", required=True,
                    help="reference proteome (Ensembl) fasta file")                    

parser.add_argument("-f", dest="fasta_file", required=True,
                    help="fasta file")                    

parser.add_argument("-log", dest="log_file", required=True,
                    help="output log file")                    

args = parser.parse_args()

print ("Reading", args.fasta_file)
fasta_entries = read_fasta(args.fasta_file)

print ("Reading", args.ref_fasta)
ref_proteins = read_fasta(args.ref_fasta)

prot_tr_ids = []
for prot in ref_proteins.values():
    trID = prot['description'].split('transcript:',1)[1].split('.',1)[0]
    prot_tr_ids.append([trID, prot['accession'].split('.',1)[0]])

tr_id_df = pd.DataFrame(data=prot_tr_ids, columns=['TranscriptID', 'ProteinID'])
prot_id_df = tr_id_df.copy()
prot_id_df.set_index('TranscriptID', inplace=True)
tr_id_df.set_index('ProteinID', inplace=True)

print ("Reading", args.input_file)
pep_df = pd.read_csv(args.input_file, sep='\t', header=0)
pep_count = len(pep_df)

# aggregated reference allele locations in protein per transcript ID
transcript_allele_locations = {}

var_db = None
if (args.var_db):
    print ("Reading", args.var_db)
    var_db = pd.read_csv(args.var_db, sep='\t', header=0)
    var_db.set_index('variantID', inplace=True)

    # add all the reference allele locations to the list
    for idx,var_row in var_db.iterrows():
        protein_canonical_allele = var_row['protein_change'].split('>',1)[0]
        protein_alt_allele = var_row['protein_change'].split('>',1)[1].split('(')[0]

        # ignore synonymous variants
        if (protein_canonical_allele == protein_alt_allele):
            continue

        if (var_row['transcriptID'] not in transcript_allele_locations):
            transcript_allele_locations[var_row['transcriptID']] = [[protein_canonical_allele], [str(var_row['chromosome']) + ':' + var_row['DNA_change'].split('>')[0]]]
        else:
            transcript_allele_locations[var_row['transcriptID']][0].append(protein_canonical_allele)
            transcript_allele_locations[var_row['transcriptID']][1].append(str(var_row['chromosome']) + ':' + var_row['DNA_change'].split('>')[0])

haplo_db = None
if (args.haplo_db):
    print ("Reading", args.haplo_db)
    haplo_db = pd.read_csv(args.haplo_db, sep='\t', header=0)
    haplo_db.set_index('HaplotypeID', inplace=True)

    # add all the reference allele locations to the list
    for idx,hap_row in haplo_db.iterrows():
        if hap_row['reading_frame'] == -1:
            continue

        # all non-synonymous alleles where the RF is known and the allele occurs after the start codon
        protein_canonical_alleles = [ var.split('>',1)[0] for var in hap_row['all_protein_changes'].split(';') if (var.split(':')[1].split('>')[0] != var.split(':')[2].split('(')[0]) and (not var.startswith('-')) ]
        DNA_alleles = [ str(hap_row['chromosome']) + ':' + hap_row['DNA_changes'].split(';')[idx].split('>')[0] for idx,var in enumerate(hap_row['all_protein_changes'].split(';')) if (var.split(':')[1].split('>')[0] != var.split(':')[2].split('(')[0]) and (not var.startswith('-')) ]
        
        if (len(protein_canonical_alleles) == 0):
            continue

        if (hap_row['TranscriptID'] not in transcript_allele_locations):
            transcript_allele_locations[hap_row['TranscriptID']] = [ protein_canonical_alleles, DNA_alleles ]
        else:
            transcript_allele_locations[hap_row['TranscriptID']][0].extend(protein_canonical_alleles)
            transcript_allele_locations[hap_row['TranscriptID']][1].extend(DNA_alleles)

# clean up the transcript allele locations
print ('Assigning reference allele locations to transcripts')
for trID in transcript_allele_locations:
    # remove duplicates
    uniq_alleles_prot = []
    uniq_alleles_DNA = []

    for idx,prot_allele in enumerate(transcript_allele_locations[trID][0]):
        if prot_allele not in uniq_alleles_prot:
            uniq_alleles_prot.append(prot_allele)
            uniq_alleles_DNA.append(transcript_allele_locations[trID][1][idx])

    positions = [ int(x.split(':')[0]) for x in uniq_alleles_prot ]
    alleles = [ x.split(':')[1] for x in uniq_alleles_prot ]
    zipped = list(zip(positions, alleles, uniq_alleles_DNA))
    zipped.sort(key=lambda x: x[0])
    try:
        positions, alleles, uniq_alleles_DNA = zip(*zipped)
    except:
        transcript_allele_locations[trID] = { 'pos': [], 'allele': [], 'DNA': [] }
        continue
    transcript_allele_locations[trID] = { 'pos': positions, 'allele': alleles, 'DNA': uniq_alleles_DNA }

# Load the annotations database
annotations_db = gffutils.FeatureDB(args.annotation_db)

log_file = open(args.log_file, 'a')
log_file.write('------------' + '[' + datetime.now().strftime('%X %x') + '] file: ' + args.input_file + ':' + '------------\n')

summary_data = []
summary_columns = ['ID', 'sequence', 'pep_type1', 'pep_type2', 'covered_changes_peptide', 'covered_changes_protein', 'covered_alleles_dna', 'matching_proteins', 'matching_transcripts', 'matching_genes', 'positions_in_proteins', 'preceding_indel_shift', 'reading_frames']

print ("Annotating peptides:")

# check if the canonical peptides contains a product of a reference allele for any of the possible variants
def check_ref_alleles(transcriptID, position, sequence, log_mismatch = False):
    if (transcriptID not in transcript_allele_locations):
        return []

    alleles = transcript_allele_locations[transcriptID]
    nearest_idx = bisect.bisect_left(alleles['pos'], position)
    result = []

    for allele_idx in range(nearest_idx, len(alleles['pos'])):
        if (alleles['pos'][allele_idx] >= position):
            if (alleles['pos'][allele_idx] >= position + len(sequence)):
                break

            pep_position = alleles['pos'][allele_idx] - position

            trimmed_allele_seq = alleles['allele'][allele_idx][0:(len(sequence) - pep_position)].replace('I', 'L')
            allele_length = len(trimmed_allele_seq)

            if (sequence[pep_position:pep_position+allele_length].replace('I', 'L') == trimmed_allele_seq):
                result.append(alleles['DNA'][allele_idx])
            else:
                if (log_mismatch):
                    log_file.write(transcriptID + ': reference allele mismatched. Expected: ' + alleles['allele'][allele_idx] + ' found: ' + sequence[:pep_position] + ' ' + sequence[pep_position:])

    return result

# check if there's a canonical peptide in the reference protein if we revert all the amino acid changes
def check_canonical_peptide(pep_seq, changes, proteinID):
    changes.sort(key=lambda x: x[0])
    ref_peptide = pep_seq
    len_diff = 0

    for ch in changes:
        ch_loc = ch[0] + len_diff
        ref_peptide = ref_peptide[:ch_loc] + ch[1] + ref_peptide[ch_loc+len(ch[2]):]

        len_diff += len(ch[2]) - len(ch[1])

    is_ref = False

    try:
        is_ref = (ref_peptide.replace('I', 'L') in ref_proteins[proteinID]['sequence'].replace('I', 'L'))
    except:
        return False

    return is_ref

def process_row(index):
    row = pep_df.iloc[index]

    fasta_accessions = re.split(r"[,;]", row['Proteins'])   # all the accessions of candidate proteins
    peptide_length = len(row['Sequence'])                   # length of the peptide sequence
    fasta_peptide_starts = [ int(pos) for pos in re.split(r"[,;]", row['Positions']) ]    # positions of the peptide within respective candidate protein 

    # crap match = any matching sequence is a contaminant
    is_contaminant = any([ 'cont' in fasta_entries[fastaID]['tag'] for fastaID in fasta_accessions ])
    if is_contaminant:
        return [row['ID'], row['Sequence'], 'contaminant', 'contaminant', '-', '-', '-', '-', '-', '-', '-', '-', '-']

    # concentrate all matching proteins (haplotype or stable protein ids)
    matching_proteins = []
    matching_transcripts = []
    reading_frames = []
    matching_protein_positions = []

    # process matching fasta entries
    for i,fastaID in enumerate(fasta_accessions):
        fasta_entry = fasta_entries[fastaID]
        matching_seq_proteins = fasta_entry['matching_proteins']    # IDs of matching protein sequences (before splitting by start and stop codons)
        prot_reading_frames = fasta_entry['reading_frames']         # RFs of matching protein sequences (known only for haplotypes)
        matching_seq_positions = fasta_entry['seq_positions']       # Positions of protein sub-sequences in the complete protein (after splitting by start and stop codons)

        for j,prot_ids in enumerate(matching_seq_proteins):
            for k,prot_id in enumerate(prot_ids):         
                if prot_id.startswith('ENSP'):
                    prot_id = prot_id.split('_', 1)[0]
                    #matching_transcripts.append(tr_id_df.loc[prot_id]['TranscriptID'])                  

                matching_proteins.append(prot_id)                   # Store protein ID
                reading_frames.append(prot_reading_frames[j][k])    # Corresponding RF

                matching_protein_positions.append(matching_seq_positions[j] + fasta_peptide_starts[i])  # Store position of the peptide in the complete protein

    # Sort everything by protein accession so that canonical matches are first (ENSP is first lexicographically)
    zipped = list(zip(matching_proteins, matching_protein_positions, reading_frames))   
    zipped.sort(key=lambda x: x[0])
    matching_proteins, matching_protein_positions, reading_frames = zip(*zipped)

    pep_type2 = ''  # type by specificity (proteoform- x protein-specific x multi-gene)

    # canonical peptide = any of the matching sequences is a canonical protein
    is_canonical = any([ 'ref' in fasta_entries[fastaID]['tag'] for fastaID in fasta_accessions ])
    if is_canonical:

        # Forget matches to variant sequences as this is a canonical peptide
        # Store the ENST ID for matching transcripts
        matching_protein_positions = [ matching_protein_positions[idx] for idx,prot_id in enumerate(matching_proteins) if prot_id.startswith('ENSP') ]
        matching_proteins = [ prot_id for prot_id in matching_proteins if prot_id.startswith('ENSP') ]
        matching_transcripts = [ tr_id_df.loc[prot_id.split('.',1)[0]]['TranscriptID'] for prot_id in matching_proteins ]

        # get gene IDs only for canonical matches (by ENSP -> ENST -> ENSG)
        matching_genes = [ annotations_db[trID.split('.', 1)[0]].attributes['gene_id'][0] for trID in matching_transcripts ]
        matching_genes = list(dict.fromkeys(matching_genes))    # remove duplicates
        if (len(matching_proteins) == 1):
            pep_type2 = 'proteoform-specific'
        elif (len(matching_genes) == 1):
            pep_type2 = 'protein-specific'
        else:
            pep_type2 = 'multi-gene'

        dna_alleles = []
        for idx,trID in enumerate(matching_transcripts):
            transcript_alleles = check_ref_alleles(trID.split('.',1)[0], matching_protein_positions[idx], row['Sequence'], True)
            if (len(transcript_alleles) > 0):
                dna_alleles.append(';'.join(transcript_alleles))
        dna_alleles = list(dict.fromkeys(dna_alleles))

        return [row['ID'], row['Sequence'], 'canonical', pep_type2, '', '', '|'.join(dna_alleles), ';'.join(matching_proteins), ';'.join(matching_transcripts), ';'.join(matching_genes), ';'.join([str(pos) for pos in matching_protein_positions]), '-', '-']

    # Here, the peptide doesn't match to any canonical (ENSP*) sequence -> annotate variation
    matching_pep_changes = []           # all unique matching changes with coordinated mapped to this peptide
    matching_protein_changes = []       # all unique matching changes in the protein
    matching_DNA_alleles = []           # corresponding unique matches to changes in the DNA
    all_preceding_indels = []               # shift of the peptide compared to the reference sequence due to preceding indels
    min_changes_found = 999999          # minimum number of found co-occurring changes (to derive the peptide type)
    has_frameshift = False              # one of the matching regions occurs after a frameshift variant
    has_canonical_alternative = False   # when all amino acid changes are reverted, the peptide matches a canonical sequence
    
    # Check if any of the peptides match to a variant outside of a haplotype 
    # check the priority haplotypes x variants
    found_variant = False

    for i,protID in enumerate(matching_proteins):
        reading_frame = int(reading_frames[i])

        # check for proteins that were added manually and don't have associated metadata
        # assume these are variants -> report protein ID as the protein change, DNA change unknown
        if protID.startswith(args.manual_prefix):
            matching_protein_changes.append(protID)
            matching_DNA_alleles.append('unknown')
            found_variant = True
            min_changes_found = min(min_changes_found, 1)

        elif ((var_db is not None) and protID.startswith(args.var_prefix)):
            variant = var_db.loc[protID]

            # there can be replicates of the same variant (i.e. the protID is not unique)
            # - in this case the .loc does not return a Series but a DataFrame with identical rows -> choose the first row of the DF if so
            if (type(variant) == type(var_db)):
                variant = variant.iloc[0]

            parent_transcript = variant['transcriptID']
            protein_start = variant['protein_prefix_length']

            if parent_transcript not in matching_transcripts:
                matching_transcripts.append(parent_transcript)

            # Peptide location within the variant sequence, accounting for the offset of the start codon (if known, first M will be at position 0)
            pep_loc = [matching_protein_positions[i] - protein_start, matching_protein_positions[i] - protein_start + peptide_length]

            local_matching_alleles_DNA = []     # corresponding alleles in the DNA
            local_matching_alleles_DNA.extend(check_ref_alleles(parent_transcript, pep_loc[0], row['Sequence']))

            protein_change = variant['protein_change']
            if ('|' in protein_change):
                protein_change = protein_change.split('|')[reading_frame]

            change_loc = int(protein_change.split('>',1)[1].split(':', 1)[0])   # location of the allele in the protein
            ref_prot_allele = protein_change.split(':', 1)[1].split('>', 1)[0].replace('I', 'L').replace('-', '')
            alt_prot_allele = protein_change.split(':', 2)[2].split('(', 1)[0].replace('I', 'L').replace('-', '')  # frameshifts are labelled with '(+fs)' at the end

            if change_loc < pep_loc[0]:
                has_frameshift = protein_change.endswith('(+fs)')

            elif ((ref_prot_allele != alt_prot_allele) or protein_change.endswith('(+fs)')) and (change_loc >= pep_loc[0]) and (change_loc < pep_loc[1]):
                # check if a frameshift happens in this peptide
                has_frameshift = has_frameshift or protein_change.endswith('(+fs)')

                # locate the alternative allele sequence in the peptide
                change_pep_loc = [ change_loc - pep_loc[0], change_loc - pep_loc[0] + len(alt_prot_allele) ]
                found_allele = row['Sequence'][change_pep_loc[0]:change_pep_loc[1]].replace('I', 'L')

                # in case the allele isn't completely covered in this peptide (i.e. there is a cleavage site in the mutated part)
                alt_prot_allele = alt_prot_allele[:len(found_allele)]

                # Sanity check: have we found the alternative allele in the peptide?
                if found_allele != alt_prot_allele:
                    log_file.write('peptide ' + row['ID'] + ' variant: ' + protID + ' expected: ' + alt_prot_allele + ' found: ' + row['Sequence'][:change_pep_loc[0]] + ' ' + row['Sequence'][change_pep_loc[0]:change_pep_loc[1]] + '\n')
                else:
                    # All looks ok -> store this change as identified
                    matching_protein_changes.append(parent_transcript + ':' + protein_change)
                    matching_pep_changes.append(str(change_pep_loc[0]) + ':' + ref_prot_allele + '>' + alt_prot_allele)

                    DNA_change_str = str(variant['chromosome']) + ':' + variant['DNA_change']
                    if DNA_change_str not in local_matching_alleles_DNA:
                        local_matching_alleles_DNA.append(DNA_change_str)

                    found_variant = True
                    min_changes_found = min(min_changes_found, 1)

            local_matching_alleles_DNA.sort(key=lambda x: int(x.split(':')[1]))
            matching_DNA_alleles.append(';'.join(local_matching_alleles_DNA))

        elif ((haplo_db is not None) and protID.startswith(args.haplo_prefix)):    

            haplotype = haplo_db.loc[protID]
            parent_transcript = haplotype['TranscriptID']
            protein_start = haplotype['protein_prefix_length']
            haplo_matching_changes = [] # the changes included in this haplotype that match the peptide

            # Remember the transcript ID for this haplotype
            if parent_transcript not in matching_transcripts:
                matching_transcripts.append(parent_transcript)

            # Peptide location within the haplotype sequence, accounting for the offset of the start codon (if known, first M will be at position 0)
            pep_loc = [matching_protein_positions[i] - protein_start, matching_protein_positions[i] - protein_start + peptide_length]
            local_matching_changes_prot = []    # hits to changes in the protein for this haplotype
            local_matching_changes_pep = []     # mapping of these changes to this peptide
            local_matching_alleles_DNA = []     # corresponding alleles in the DNA
            preceding_indels = 0                # shift of the peptide compared to the reference sequence due to preceding indels

            all_protein_changes = haplotype['all_protein_changes'].split(';')
            all_DNA_changes = haplotype['DNA_changes'].split(';')

            # if there are multiple possible reading frames, check only the changes in the matching reading frame
            if ('|' in haplotype['all_protein_changes']):
                all_protein_changes = [ ch.split('|')[reading_frame] for ch in all_protein_changes ]        

            for j,ch in enumerate(all_protein_changes):
                change_loc = int(ch.split('>', 1)[1].split(':',1)[0])   # location of the alternative allele in the protein (can be shifted by preceding indels in the haplotype)
                ref_prot_allele = ch.split(':', 1)[1].split('>', 1)[0].replace('I', 'L').replace('-', '')
                alt_prot_allele = ch.split(':', 2)[2].split('(', 1)[0].replace('I', 'L').replace('-', '')  # frameshifts are labelled with '(+fs)' at the end

                if change_loc < pep_loc[0]:
                    has_frameshift = has_frameshift or ch.endswith('(+fs)')
                    preceding_indels += len(alt_prot_allele) - len(ref_prot_allele)
                
                # is it not a synonymous mutation, and does it happen in this peptide?
                elif ((ref_prot_allele != alt_prot_allele) or ch.endswith('(+fs)')) and (change_loc >= pep_loc[0]) and (change_loc < pep_loc[1]):
                    # check if a frameshift happens in this peptide
                    has_frameshift = has_frameshift or ch.endswith('(+fs)')

                    # locate the alternative allele sequence in the peptide
                    change_pep_loc = [ change_loc - pep_loc[0], change_loc - pep_loc[0] + len(alt_prot_allele) ]
                    found_allele = row['Sequence'][change_pep_loc[0]:change_pep_loc[1]].replace('I', 'L')

                    # in case the allele isn't completely covered in this peptide (i.e. there is a cleavage site in the mutated part)
                    alt_prot_allele = alt_prot_allele[:len(found_allele)]

                    # Sanity check: have we found the alternative allele in the peptide?
                    if found_allele != alt_prot_allele:
                        log_file.write('PSM ' + row['ID'] + ' haplotype:' + protID + ' expected:' + alt_prot_allele + ' found: ' + row['Sequence'][:change_pep_loc[0]] + ' ' + row['Sequence'][change_pep_loc[0]:])
                    else:
                        # All looks ok -> store this change as identified
                        haplo_matching_changes.append([change_pep_loc[0], ref_prot_allele, alt_prot_allele])
                        local_matching_changes_prot.append(ch)
                        local_matching_changes_pep.append(str(change_pep_loc[0]) + ':' + ref_prot_allele + '>' + alt_prot_allele)
                        chromosome = protID.split('chr',1)[1].split('_',1)[0]
                        local_matching_alleles_DNA.append(chromosome + ':' + all_DNA_changes[j])

            haplo_has_canonical_alternative = check_canonical_peptide(row['Sequence'], haplo_matching_changes, prot_id_df.loc[parent_transcript]['ProteinID'])
            has_canonical_alternative = has_canonical_alternative or haplo_has_canonical_alternative

            # Update the minimal number of co-occurring changes found
            # Note: if we find a peptide that matches to a canonical region with variation and to a completely novel region with no variation, it is still a variant peptide
            if (haplo_has_canonical_alternative):
                min_changes_found = min(min_changes_found, len(local_matching_alleles_DNA))

            # get the reference alleles covered
            local_matching_alleles_DNA.extend(check_ref_alleles(parent_transcript, pep_loc[0] - preceding_indels, row['Sequence']))
            local_matching_alleles_DNA.sort(key=lambda x: int(x.split(':')[1]))
            DNA_alleles_str = ';'.join(local_matching_alleles_DNA)

            all_preceding_indels.append(preceding_indels)

            # if two variants affect the same codon, the change will be shown twice -> remove the duplicate
            local_matching_changes_pep = list(dict.fromkeys(local_matching_changes_pep))
            local_matching_changes_prot = list(dict.fromkeys(local_matching_changes_prot))

            if (len(local_matching_changes_prot) > 0):
                prot_changes_str = parent_transcript + ':' + ';'.join(local_matching_changes_prot)
                matching_pep_changes.append(';'.join(local_matching_changes_pep))
            else:
                prot_changes_str = ''

            if prot_changes_str not in matching_protein_changes:
                matching_protein_changes.append(prot_changes_str)
            if DNA_alleles_str not in matching_DNA_alleles:
                matching_DNA_alleles.append(DNA_alleles_str)

        else: # this should not happen!
            raise ValueError('Protein ID in unknown format:', protID)

    # Get  gene IDs from matching transcripts (ENST -> ENSG)
    matching_genes = [ annotations_db[trID].attributes['gene_id'][0] for trID in matching_transcripts ]
    matching_genes = list(dict.fromkeys(matching_genes))    # remove duplicate gene IDs
    matching_DNA_alleles = list(dict.fromkeys(matching_DNA_alleles)) # remove duplicate combinations of alleles
    matching_pep_changes = list(dict.fromkeys(matching_pep_changes))

    # Check if we have found any alternative alleles. If not, set the alleles counter to 0.
    has_alt_allele = any([ ('>' in alleles_str) for alleles_str in matching_DNA_alleles ])

    if (len(matching_proteins) == 1):
        pep_type2 = 'proteoform-specific'
    elif (len(matching_genes) == 1):
        pep_type2 = 'protein-specific'
    else:
        pep_type2 = 'multi-gene'

    if found_variant:
        return [row['ID'], row['Sequence'], 'single-variant(ProVar)', pep_type2, '|'.join(matching_pep_changes), '|'.join(matching_protein_changes), '|'.join(matching_DNA_alleles), ';'.join(matching_proteins), ';'.join(matching_transcripts), ';'.join(matching_genes), ';'.join([str(pos) for pos in matching_protein_positions]), ';'.join([str(x) for x in all_preceding_indels]) if (len(all_preceding_indels) > 0) else '-',  ';'.join(reading_frames)]
    if (min_changes_found > 1) and has_canonical_alternative:
        return [row['ID'], row['Sequence'], 'multi-variant', pep_type2, '|'.join(matching_pep_changes), '|'.join(matching_protein_changes), '|'.join(matching_DNA_alleles), ';'.join(matching_proteins), ';'.join(matching_transcripts), ';'.join(matching_genes), ';'.join([str(pos) for pos in matching_protein_positions]), ';'.join([str(x) for x in all_preceding_indels]) if (len(all_preceding_indels) > 0) else '-', ';'.join(reading_frames)]
    elif has_alt_allele and has_canonical_alternative:
        return [row['ID'], row['Sequence'], 'single-variant', pep_type2, '|'.join(matching_pep_changes), '|'.join(matching_protein_changes), '|'.join(matching_DNA_alleles), ';'.join(matching_proteins), ';'.join(matching_transcripts), ';'.join(matching_genes), ';'.join([str(pos) for pos in matching_protein_positions]), ';'.join([str(x) for x in all_preceding_indels]) if (len(all_preceding_indels) > 0) else '-', ';'.join(reading_frames)]
    elif has_alt_allele and not has_canonical_alternative and not has_frameshift:
        return [row['ID'], row['Sequence'], 'variant-no-ref', pep_type2, '|'.join(matching_pep_changes), '|'.join(matching_protein_changes), '|'.join(matching_DNA_alleles), ';'.join(matching_proteins), ';'.join(matching_transcripts), ';'.join(matching_genes), ';'.join([str(pos) for pos in matching_protein_positions]), ';'.join([str(x) for x in all_preceding_indels]) if (len(all_preceding_indels) > 0) else '-', ';'.join(reading_frames)]
    elif has_frameshift:
        # Sequence contains a result of a frameshift, the mutation occurs either in the peptide, or upstream in the protein before it
        return [row['ID'], row['Sequence'], 'frameshift', pep_type2, '|'.join(matching_pep_changes), '|'.join(matching_protein_changes), '|'.join(matching_DNA_alleles), ';'.join(matching_proteins), ';'.join(matching_transcripts), ';'.join(matching_genes), ';'.join([str(pos) for pos in matching_protein_positions]), ';'.join([str(x) for x in all_preceding_indels]) if (len(all_preceding_indels) > 0) else '-', ';'.join(reading_frames)]
    else:   
        # looks like a translation of a canonical CDS sequence that doesn't have a canonical protein in Ensembl (e.g., alternative reading frame)
        return [row['ID'], row['Sequence'], 'canonical-no-ref', pep_type2, '|'.join(matching_pep_changes), '|'.join(matching_protein_changes), '|'.join(matching_DNA_alleles), ';'.join(matching_proteins), ';'.join(matching_transcripts), ';'.join(matching_genes), ';'.join([str(pos) for pos in matching_protein_positions]), ';'.join([str(x) for x in all_preceding_indels]) if (len(all_preceding_indels) > 0) else '-', ';'.join(reading_frames)]

# store results
with Pool(args.threads) as p:
    summary_data = list(tqdm(p.imap_unordered(process_row, range(0, pep_count)), total=pep_count))
    #summary_data = list(map(process_row, range(0, pep_count)))
    p.close()
    p.join()
    summary_df = pd.DataFrame(columns=summary_columns, data=summary_data)
    summary_df.to_csv(args.output_file, header=True, index=False, sep='\t')

log_file.close()
