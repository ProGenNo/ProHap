import argparse
import pandas as pd
import re
from datetime import datetime
from tqdm import tqdm
from multiprocessing import Pool
from modules.common import read_fasta

parser = argparse.ArgumentParser(
	description='Reads the PSM report file, creates a new file with information about covered SNPs for each PSM and protein')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input CSV file")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

parser.add_argument("-id_col", dest="id_column", required=False,
                    help="ID columns name, default: PSMId", default='PSMId')

parser.add_argument("-hap", dest="haplo_db", required=False,
                    help="haplotypes tab-separated file (optional)", default=None)

parser.add_argument("-var", dest="var_db", required=False,
                    help="variants tab-separated file (optional)", default=None)

parser.add_argument("-hap_prefix", dest="haplo_prefix", required=False,
                    help="prefix for haplotype protein ID (default: 'haplo_')", default='haplo_')

parser.add_argument("-var_prefix", dest="var_prefix", required=False,
                    help="prefix for variant protein ID (default: 'var_')", default='var_')

parser.add_argument("-manual_prefix", dest="manual_prefix", required=False,
                    help="prefix for manually added protein sequences (default: 'man_')", default='man_')

parser.add_argument("-tr_id", dest="transcript_ids", required=True,
                    help="csv file mapping protein IDs to transcript IDs")

parser.add_argument("-g_id", dest="gene_ids", required=True,
                    help="csv file mapping transcript IDs to gene IDs")

parser.add_argument("-t", dest="threads", type=int, required=True,
                    help="maximum number of threads")

parser.add_argument("-f", dest="fasta_file", required=True,
                    help="fasta file")                    

parser.add_argument("-log", dest="log_file", required=True,
                    help="output log file")                    

args = parser.parse_args()

print ("Reading", args.fasta_file)
fasta_entries = read_fasta(args.fasta_file)

print ("Reading", args.input_file)
psm_df = pd.read_csv(args.input_file, sep='\t', header=0)
psm_count = len(psm_df)

var_db = None
if (args.var_db):
    print ("Reading", args.var_db)
    var_db = pd.read_csv(args.var_db, sep='\t', header=0)
    var_db.set_index('variantID', inplace=True)

haplo_db = None
if (args.haplo_db):
    print ("Reading", args.haplo_db)
    haplo_db = pd.read_csv(args.haplo_db, sep='\t', header=0)
    haplo_db.set_index('HaplotypeID', inplace=True)

print ("Reading", args.transcript_ids)
tr_id_df = pd.read_csv(args.transcript_ids, header=0)
tr_id_df.set_index('ProteinID', inplace=True)

print ("Reading", args.gene_ids)
gene_id_df = pd.read_csv(args.gene_ids, header=0)
gene_id_df.set_index('TranscriptID', inplace=True)

log_file = open(args.log_file, 'a')
log_file.write('------------' + '[' + datetime.now().strftime('%X %x') + '] file: ' + args.input_file + ':' + '------------\n')

summary_data = []
summary_columns = [args.id_column, 'psm_type1', 'psm_type2', 'covered_changes_protein', 'covered_changes_dna', 'matching_proteins', 'matching_transcripts', 'matching_genes', 'positions_in_proteins', 'reading_frames']

print ("Annotating PSMs:")

def process_row(index):
    row = psm_df.iloc[index]

    fasta_accessions = [ acc.split('.')[0] for acc in re.split(r"[,;]", row['Proteins']) ]  # all the accessions of candidate proteins
    peptide_length = len(row['Sequence'])                   # length of the peptide sequence
    fasta_peptide_starts = [ int(pos) for pos in re.split(r"[,;]", row['Positions']) ]    # positions of the peptide within respective candidate protein 

    # decoy match = all matching sequences are decoys
    is_decoy = all([ fastaID.endswith('REVERSED') for fastaID in fasta_accessions ])
    if is_decoy:
        return [row[args.id_column], 'decoy', 'decoy', '', '', '', '', '', '', '']

    # filter out matching decoy elements in fasta (in case there is an overlap between a decoy and target sequence), and remember only corresponding peptide positions
    fasta_peptide_starts = [ fasta_peptide_starts[i] for i,fastaID in enumerate(fasta_accessions) if not fastaID.endswith('REVERSED') ]
    fasta_accessions = [ fastaID for fastaID in fasta_accessions if not fastaID.endswith('REVERSED') ]

    # crap match = any matching sequence is a contaminant
    is_contaminant = any([ 'cont' in fasta_entries[fastaID]['tag'] for fastaID in fasta_accessions ])
    if is_contaminant:
        return [row[args.id_column], 'contaminant', 'contaminant', '', '', '', '', '', '', '']

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
                    prot_id = prot_id.split('_', 1)[0].strip()
                    matching_transcripts.append(tr_id_df.loc[prot_id]['TranscriptID'])                  # Transcript ID only for canonical sequences (where ENSP is known)
                try:
                    matching_proteins.append(prot_id)                   # Store protein ID
                    reading_frames.append(prot_reading_frames[j][k])    # Corresponding RF
                except:
                    print (prot_id)
                    print (reading_frames)

                matching_protein_positions.append(matching_seq_positions[j] + fasta_peptide_starts[i])  # Store position of the peptide in the complete protein

    # Sort everything by protein accession so that canonical matches are first (ENSP is first lexicographically)
    zipped = list(zip(matching_proteins, matching_protein_positions, reading_frames))   
    zipped.sort(key=lambda x: x[0])
    matching_proteins, matching_protein_positions, reading_frames = zip(*zipped)

    psm_type2 = ''  # type by specificity (proteoform- x protein-specific x multi-gene)

    # canonical PSM = any of the matching sequences is a canonical protein
    is_canonical = any([ 'ref' in fasta_entries[fastaID]['tag'] for fastaID in fasta_accessions ])
    if is_canonical:

        # get gene IDs only for canonical matches (by ENSP -> ENST -> ENSG)
        matching_genes = [ gene_id_df.loc[trID.split('.', 1)[0]]['GeneID'] for trID in matching_transcripts ]
        matching_genes = list(dict.fromkeys(matching_genes))    # remove duplicates
        if (len(matching_proteins) == 1):
            psm_type2 = 'proteoform-specific'
        elif (len(matching_genes) == 1):
            psm_type2 = 'protein-specific'
        else:
            psm_type2 = 'multi-gene'

        return [row[args.id_column], 'canonical', psm_type2, '', '', ';'.join(matching_proteins), ';'.join(matching_transcripts), ';'.join(matching_genes),';'.join([str(pos) for pos in matching_protein_positions]), ';'.join(reading_frames)]

    # Here, the peptide doesn't match to any canonical (ENSP*) sequence -> annotate variation
    matching_protein_changes = []   # all unique matching changes in the protein
    matching_DNA_changes = []       # corresponding unique matches to changes in the DNA
    min_changes_found = 999999      # minimum number of found co-occurring changes (to derive the peptide type)
    has_frameshift = False          # one of the matching regions occurs after a frameshift variant
    
    # Check if any of the PSMs match to a variant outside of a haplotype 
    # TODO: check the priority haplotypes x variants
    found_variant = False

    for i,protID in enumerate(matching_proteins):
        reading_frame = int(reading_frames[i])

        # check for proteins that were added manually and don't have associated metadata
        # assume these are variants -> report protein ID as the protein change, DNA change unknown
        if protID.startswith(args.manual_prefix):
            matching_protein_changes.append(protID)
            matching_DNA_changes.append('unknown')
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

            protein_change = variant['protein_change']
            if ('|' in protein_change):
                protein_change = protein_change.split('|')[reading_frame]

            change_loc = int(protein_change.split('>',1)[1].split(':', 1)[0])   # location of the allele in the protein
            ref_prot_allele = protein_change.split(':', 1)[1].split('>', 1)[0].replace('I', 'L').replace('-', '')
            alt_prot_allele = protein_change.split(':', 2)[2].split('(', 1)[0].replace('I', 'L').replace('-', '')  # frameshifts are labelled with '(+fs)' at the end

            if (ref_prot_allele != alt_prot_allele) and (change_loc >= pep_loc[0]) and (change_loc < pep_loc[1]):

                # locate the alternative allele sequence in the peptide
                change_pep_loc = [ change_loc - pep_loc[0], change_loc - pep_loc[0] + len(alt_prot_allele) ]
                found_allele = row['Sequence'][change_pep_loc[0]:change_pep_loc[1]].replace('I', 'L')

                # in case the allele isn't completely covered in this peptide (i.e. there is a cleavage site in the mutated part)
                alt_prot_allele = alt_prot_allele[:len(found_allele)]

                # Sanity check: have we found the alternative allele in the peptide?
                if found_allele != alt_prot_allele:
                    log_file.write('PSM ' + row[args.id_column] + ' variant:' +  protID + ' expected:' + alt_prot_allele + ' found:' + row['Sequence'][:change_pep_loc[0]] + ' ' + row['Sequence'][change_pep_loc[0]:change_pep_loc[1]] + '\n')
                else:
                # All looks ok -> store this change as identified
                    matching_protein_changes.append(parent_transcript + ':' + protein_change)

                    DNA_change_str = str(variant['chromosome']) + ':' + variant['DNA_change']
                    if DNA_change_str not in matching_DNA_changes:
                        matching_DNA_changes.append(DNA_change_str)

                    found_variant = True
                    min_changes_found = min(min_changes_found, 1)

            else:
                # check if frameshift happened before
                if change_loc < pep_loc[0]:
                    has_frameshift = protein_change.endswith('+fs')

                # remember that that this is a translation of a canonical cDNA    
                min_changes_found = 0

        elif ((haplo_db is not None) and protID.startswith(args.haplo_prefix)):    

            haplotype = haplo_db.loc[protID]
            parent_transcript = haplotype['TranscriptID']
            protein_start = haplotype['protein_prefix_length']

            # Remember the transcript ID for this haplotype
            if parent_transcript not in matching_transcripts:
                matching_transcripts.append(parent_transcript)

            # Peptide location within the haplotype sequence, accounting for the offset of the start codon (if known, first M will be at position 0)
            pep_loc = [matching_protein_positions[i] - protein_start, matching_protein_positions[i] - protein_start + peptide_length]
            local_matching_changes_prot = []    # hits to changes in the protein for this haplotype
            local_matching_changes_DNA = []     # corresponding changes in the DNA

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
                    has_frameshift = has_frameshift or ch.endswith('+fs')
                
                # is it not a synonymous mutation, and does it happen in this peptide?
                elif (ref_prot_allele != alt_prot_allele) and (change_loc >= pep_loc[0]) and (change_loc < pep_loc[1]):

                    # locate the alternative allele sequence in the peptide
                    change_pep_loc = [ change_loc - pep_loc[0], change_loc - pep_loc[0] + len(alt_prot_allele) ]
                    found_allele = row['Sequence'][change_pep_loc[0]:change_pep_loc[1]].replace('I', 'L')

                    # in case the allele isn't completely covered in this peptide (i.e. there is a cleavage site in the mutated part)
                    alt_prot_allele = alt_prot_allele[:len(found_allele)]

                    # Sanity check: have we found the alternative allele in the peptide?
                    if found_allele != alt_prot_allele:
                        log_file.write('PSM ' + row[args.id_column] + ' haplotype:' + protID + ' expected:' + alt_prot_allele + ' found:' + row['Sequence'][:change_pep_loc[0]] + ' ' + row['Sequence'][change_pep_loc[0]:])
                    else:
                        # All looks ok -> store this change as identified
                        local_matching_changes_prot.append(ch)
                        chromosome = protID.split('chr',1)[1].split('_',1)[0]
                        local_matching_changes_DNA.append(chromosome + ':' + all_DNA_changes[j])

            if (len(local_matching_changes_prot) > 0):
                prot_changes_str = parent_transcript + ':' + ';'.join(local_matching_changes_prot)
                DNA_changes_str = ';'.join(local_matching_changes_DNA)
            else:
                prot_changes_str = ''
                DNA_changes_str = ''

            # Update the minimal number of co-occurring changes found
            min_changes_found = min(min_changes_found, len(local_matching_changes_DNA))

            if prot_changes_str not in matching_protein_changes:
                matching_protein_changes.append(prot_changes_str)
            if DNA_changes_str not in matching_DNA_changes:
                matching_DNA_changes.append(DNA_changes_str)

        else: # this should not happen!
            raise ValueError('Protein ID in unknown format:', protID)

    # Get  gene IDs from matching transcripts (ENST -> ENSG)
    matching_genes = [ gene_id_df.loc[trID]['GeneID'] for trID in matching_transcripts ]
    matching_genes = list(dict.fromkeys(matching_genes))    # remove duplicates
    if (len(matching_proteins) == 1):
        psm_type2 = 'proteoform-specific'
    elif (len(matching_genes) == 1):
        psm_type2 = 'protein-specific'
    else:
        psm_type2 = 'multi-gene'

    if found_variant:
        return [row[args.id_column], 'vcf-variant', psm_type2, '|'.join(matching_protein_changes), '|'.join(matching_DNA_changes), ';'.join(matching_proteins), ';'.join(matching_transcripts), ';'.join(matching_genes), ';'.join([str(pos) for pos in matching_protein_positions]), ';'.join(reading_frames)]
    if (min_changes_found > 1):
        return [row[args.id_column], 'multi-variant', psm_type2, '|'.join(matching_protein_changes), '|'.join(matching_DNA_changes), ';'.join(matching_proteins), ';'.join(matching_transcripts), ';'.join(matching_genes), ';'.join([str(pos) for pos in matching_protein_positions]), ';'.join(reading_frames)]
    elif (min_changes_found > 0):
        return [row[args.id_column], 'single-variant', psm_type2, '|'.join(matching_protein_changes), '|'.join(matching_DNA_changes), ';'.join(matching_proteins), ';'.join(matching_transcripts), ';'.join(matching_genes), ';'.join([str(pos) for pos in matching_protein_positions]), ';'.join(reading_frames)]
    elif has_frameshift:
        # in some of the matching proteins, no changes were mapped to this peptide, but it occurs after a frameshift, so no canonical protein was matched
        return [row[args.id_column], 'canonical-frameshift', psm_type2, '|'.join(matching_protein_changes), '|'.join(matching_DNA_changes), ';'.join(matching_proteins), ';'.join(matching_transcripts), ';'.join(matching_genes), ';'.join([str(pos) for pos in matching_protein_positions]), ';'.join(reading_frames)]
    else:   
        # looks like a translation of a canonical CDS sequence that doesn't have a canonical protein in Ensembl (e.g., alternative reading frame)
        return [row[args.id_column], 'canonical-no-ref', psm_type2, '|'.join(matching_protein_changes), '|'.join(matching_DNA_changes), ';'.join(matching_proteins), ';'.join(matching_transcripts), ';'.join(matching_genes), ';'.join([str(pos) for pos in matching_protein_positions]), ';'.join(reading_frames)]

# store results
with Pool(args.threads) as p:
    summary_data = list(tqdm(p.imap_unordered(process_row, range(0, psm_count)), total=psm_count))
#    summary_data = list(map(process_row, range(0, psm_count)))
    p.close()
    p.join()
    summary_df = pd.DataFrame(columns=summary_columns, data=summary_data)
    summary_df.to_csv(args.output_file, header=True, index=False)

log_file.close()
