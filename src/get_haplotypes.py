import pandas as pd
import bisect

result_columns = ['TranscriptID', 'Changes', 'AlleleFrequencies', 'RemovedChanges', 'VCF_IDs', 'Count', 'Frequency', 'Samples']
result_data = []

# group mutations to see which ones conflict
def cluster_conflicting_mutations(changes):
    eventQ = [ { 'loc': ch['change'][0], 'type': 's', 'id': ch['id'] } for ch in changes ]
    eventQ.extend([ { 'loc': ch['change'][0] + len(ch['change'][1]), 'type': 'e', 'id': ch['id'] } for ch in changes ])

    eventQ.sort(key=lambda x: x['loc'])

    active_ids = []     # mutations overlapping current position
    id_groups = []      # all the groups of mutations already passed
    current_group = []  # list for aggregating currently overlapping mutations

    for evt in eventQ:
        if evt['type'] == 's':
            active_ids.append(evt['id'])
            current_group.append(evt['id'])

        elif evt['type'] == 'e':
            active_ids.remove(evt['id'])

            if (len(active_ids) == 0):
                id_groups.append(current_group)
                current_group = []

    return id_groups

# check the list of mutations for any potential conflicts (multiple mutations affecting the same locus)
# if there are multiple mutations conflicting, gradually remove the one with lowest AF until no conflicts
def remove_conflicting_mutations(changes, AFs):    
    changes_enum = [{ 'change': ch, 'id': i} for i,ch in enumerate(changes)]
    id_groups = cluster_conflicting_mutations(changes_enum)

    result_kept = []
    result_removed = []

    while (len(id_groups) > 0):
        tmp_groups = []

        for group in id_groups:
            if (len(group) == 1):               # no conflifting mutations here
                result_kept.append(group[0])
            else:
                group_sorted = sorted(group, key=lambda i: float(AFs[i]))
                result_removed.append(group_sorted.pop(0)) # remove the mutation with lowest AF
                group_changes = [ changes_enum[i] for i in group_sorted ]

                tmp_groups.extend(cluster_conflicting_mutations(group_changes))

        id_groups = tmp_groups

    '''
    for group in id_groups:
        while (len(group) > 1)


        if (len(group) == 1):               # no conflifting mutations here
            result_kept.append(group[0])
        elif (len(group) > 1):              # conflicting mutations -> sort according to AF and pick the highest one, remove the rest
            group_sorted = sorted(group, key=lambda i: float(AFs[i]))
            result_kept.append(group_sorted[0])
            result_removed.extend(group_sorted[1:])
    '''
    return result_kept, result_removed

# Creates a list of observed haplotypes from VCF files (individual file for each transcript, with phased genotypes). 
# Returns a dataframe, haplotypes described by DNA location, reference and alternative allele.
def get_gene_haplotypes(all_transcripts, vcf_dfs):
    # check haplotypes for every transcript in the DB
    for transcript_idx,transcript in enumerate(all_transcripts):
        transcriptID = transcript.id
        
        # load the according VCF file to Pandas (skip the comments at the beginning)
        vcf_df = vcf_dfs[transcriptID]

        # IDs of phased genotype columns
        indiv_ids = vcf_df.columns.values[9:]
        indiv_count = len(indiv_ids)

        # no variation in this transcript -> store reference haplotype only
        if (len(vcf_df) == 0):
            result_data.append([transcriptID, 'REF', '', '', 'REF', indiv_count * 2, 1.0, 'all'])
            continue

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

    return result_df

    #result_df.to_csv(args.output_file, sep='\t', header=True, index=False)
    
# checksum
# print(sum(result_df['Frequency'].to_list()))
