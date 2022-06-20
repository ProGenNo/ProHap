import pandas as pd
import argparse
import gffutils
import bisect

samples_file = "igsr_samples.tsv"

parser = argparse.ArgumentParser(
    description="Creates a list of observed haplotypes from VCF files (individual file for each transcript, with phased genotypes). Exported as a CSV, haplotypes described by DNA location, reference and alternative allele.")

parser.add_argument("-db", dest="annotation_db", required=True,
                    help="DB file created by gffutils from GTF")

parser.add_argument("-d", dest="input_dir", required=True,
                    help="input directory")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output DB")

args = parser.parse_args()

#samples_df = pd.read_csv(samples_file, sep='\t', index_col=0)

result_columns = ['TranscriptID', 'Changes', 'AlleleFrequencies', 'RemovedChanges', 'VCF_IDs', 'Count', 'Frequency', 'Samples']
result_data = []

# Load the annotations database
annotations_db = gffutils.FeatureDB(args.annotation_db)

all_transcripts = [ transcript for transcript in annotations_db.features_of_type('transcript', order_by='start') ]

# check the list of mutations for any potential conflicts (multiple mutations affecting the same locus)
# of there are multiple mutations conflicting, keep the one with highest AF
def remove_conflicting_mutations(changes, AFs):    
    eventQ = [ { 'loc': ch[0], 'type': 's', 'id': i } for i,ch in enumerate(changes) ]
    eventQ.extend([ { 'loc': ch[0] + len(ch[1]), 'type': 'e', 'id': i } for i,ch in enumerate(changes) ])

    eventQ.sort(key=lambda x: x['loc'])

    active_ids = []     # mutations overlapping current position
    id_groups = []      # all the groups of mutations already passed
    current_group = []  # list for aggregating currently overlapping mutations

    result_kept = []
    result_removed = []

    for evt in eventQ:
        if evt['type'] == 's':
            active_ids.append(evt['id'])
            current_group.append(evt['id'])

        elif evt['type'] == 'e':
            active_ids.remove(evt['id'])

            if (len(active_ids) == 0):
                id_groups.append(current_group)
                current_group = []

    for group in id_groups:
        if (len(group) == 1):               # no conflifting mutations here
            result_kept.append(group[0])
        elif (len(group) > 1):              # conflicting mutations -> sort according to AF and pick the highest one, remove the rest
            group_sorted = sorted(group, key=lambda i: -float(AFs[i]))
            result_kept.append(group_sorted[0])
            result_removed.extend(group_sorted[1:])

    return result_kept, result_removed

# check haplotypes for every transcript in the DB
for transcript_idx,transcript in enumerate(all_transcripts):
    transcriptID = transcript.id
    
    # load the according VCF file to Pandas (skip the comments at the beginning)
    vcf_df = pd.read_csv(args.input_dir + '/' + transcriptID + '.vcf', sep='\t', header=19)
    if (len(vcf_df) == 0):
        continue

    # IDs of phased genotype columns
    indiv_ids = vcf_df.columns.values[9:]
    indiv_count = len(indiv_ids)

    haplo_combinations = []
    haplo_samples = []

    # check the combination for every individual
    for indiv in indiv_ids:

        # store indices of rows for which the alternative allele has been found -> create a temporary string ID of the haplotype
        vals = vcf_df[indiv].to_list()

        # sanity check - correct separator between paternal / maternal chromosome
        err_rows = ','.join([ str(i) for i,elem in enumerate(vals) if '|' not in elem ])
        if (len(err_rows) > 1):
            print('Incorrect formatting!', 'indivudial:', indiv, 'rows:', err_rows, 'transcript:', transcriptID)

        hap1 = ','.join([ str(i) for i,elem in enumerate(vals) if elem.startswith('1|') ])
        hap2 = ','.join([ str(i) for i,elem in enumerate(vals) if elem.endswith('|1') ])

        # no alternative alleles -> reference haplotype
        if hap1 == '':
            hap1 = 'REF'
        if hap2 == '':
            hap2 = 'REF'

        # find if this haplotype has been identified before
        # if so, add this individual to the list, otherwise add new haoplotype
        nearest_idx = bisect.bisect_left(haplo_combinations, hap1)
        if (len(haplo_combinations) > nearest_idx and haplo_combinations[nearest_idx] == hap1):
            haplo_samples[nearest_idx].append(indiv + ':1')
        else:
            haplo_combinations.insert(nearest_idx, hap1)
            haplo_samples.insert(nearest_idx, [indiv + ':1'])

        nearest_idx = bisect.bisect_left(haplo_combinations, hap2)
        if (len(haplo_combinations) > nearest_idx and haplo_combinations[nearest_idx] == hap2):
            haplo_samples[nearest_idx].append(indiv + ':2')
        else:
            haplo_combinations.insert(nearest_idx, hap2)
            haplo_samples.insert(nearest_idx, [indiv + ':2'])    

    # once all individuals in this VCF have been processed -> summarize observed haplotypes, compute worldwide frequencies
    for i,combination in enumerate(haplo_combinations):
        if combination == 'REF':
            changes_str = 'REF'
            AFs_str = ""
            combination = ""
            removed_str = ""

        else:
            indexes = [ int(idx) for idx in combination.split(',') ]
            changes = []    # changes in the POS:REF>ALT format
            changelist = [] # changes with the POS, REF and ALT fields sepatared
            AFs = []        # allele frequencies
            for idx in indexes:
                row = vcf_df.iloc[idx]

                changelist.append([row['POS'], row['REF'], row['ALT']])
                changes.append(str(row['POS']) + ':' + row['REF'] + '>' + row['ALT'])
                if 'AF' in row['INFO']:
                    AFs.append(row['INFO'].split('AF=')[1].split(';')[0])
                else:
                    AFs.append('-1')

            # check for conflicting mutations!
            kept, removed = remove_conflicting_mutations(changelist, AFs)
            removedChanges = [ changes[i] for i in removed ]
            changes = [ changes[i] for i in kept ]
            AFs = [ AFs[i] for i in kept ]

            removed_str = ';'.join(removedChanges)
            changes_str = ';'.join(changes)
            AFs_str = ';'.join(AFs)

        result_data.append([transcriptID, changes_str, AFs_str, removed_str, haplo_combinations[i], len(haplo_samples[i]), len(haplo_samples[i]) / (indiv_count * 2), ';'.join(haplo_samples[i])])

    # print(transcriptID + ': ' + str(transcript_idx) + ' / ' + str(len(all_transcripts)), end='\r')

result_df = pd.DataFrame(columns=result_columns, data=result_data)
result_df.sort_values(by=['TranscriptID', 'Frequency'], ascending=[True, False], inplace=True)
result_df.to_csv(args.output_file, sep='\t', header=True, index=False)
    
# checksum
# print(sum(result_df['Frequency'].to_list()))
