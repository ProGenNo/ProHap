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

result_columns = ['TranscriptID', 'Changes', 'Count', 'Frequency', 'Samples']
result_data = []

# Load the annotations database
annotations_db = gffutils.FeatureDB(args.annotation_db)

all_transcripts = [ transcript for transcript in annotations_db.features_of_type('transcript', order_by='start') ]

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
        hap1 = ','.join([ str(i) for i,elem in enumerate(vals) if elem.startswith('1') ])
        hap2 = ','.join([ str(i) for i,elem in enumerate(vals) if elem.endswith('1') ])

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
        else:
            indexes = [ int(idx) for idx in combination.split(',') ]
            changes = []
            for idx in indexes:
                row = vcf_df.iloc[idx]

                changes.append(str(row['POS']) + ':' + row['REF'] + '>' + row['ALT'])
            changes_str = ';'.join(changes)

        result_data.append([transcriptID, changes_str, len(haplo_samples[i]), len(haplo_samples[i]) / (indiv_count * 2), ';'.join(haplo_samples[i])])

    # print(transcriptID + ': ' + str(transcript_idx) + ' / ' + str(len(all_transcripts)), end='\r')

result_df = pd.DataFrame(columns=result_columns, data=result_data)
result_df.sort_values(by=['TranscriptID', 'Frequency'], ascending=[True, False], inplace=True)
result_df.to_csv(args.output_file, sep='\t', header=True, index=False)
    
# checksum
# print(sum(result_df['Frequency'].to_list()))
