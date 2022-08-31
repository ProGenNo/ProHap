import bisect
from common import KeyWrapper
from io import StringIO
import pandas as pd
import re

# Process a VCF file, select rows that intersect exons of given transcripts. Results are written as TSV files in to a temporary folder. Returns a list of column names in the VCF.
# input: 
# all_transcripts: list of GTF transcript features, ordered by start position
# vcf_file: file handle for reading the VCF
# annotations_db: FeatureDB of the GTF file
# min_af: threshold allele frequency (float)
def parse_vcf(all_transcripts, vcf_file, annotations_db, min_af, tmp_dir):

    # read the header of the VCF - keep only the last line of the header
    VCF_header = ""

    vcf_linecount = 1
    line = vcf_file.readline()

    while (line != "" and line.startswith('#')):
        VCF_header = line[1:]
        vcf_linecount += 1
        line = vcf_file.readline()

    # browse the chromosome in a sweep-line approach - assumes that the VCF file is sorted!
    # keep a list of transcripts that intersect the current position of the sweep line -> assign the VCF line to all of these transcripts

    # TODO: get the coordinates within the transcript already here?

    transcript_queue = []               # queue of transcript objects inc. the exons, sorted by end position, each element aggregates the VCF file contents
    current_pos = int(line.split()[1])  # position of the current VCF entry 
    #result_dfs = {}                     # a list of dataframes with VCF entries for each transcript, accessed by the stable transcript id

    last_transcript = None

    # iterate through all the transcripts - add the first one to the queue (and all others starting at the same location), and of each other, check if there is a gap that can be filled in by VCF entries
    # i.e., process all the VCF entries that lay before the current transcript -> update the queue
    for current_transcript in all_transcripts:

        if (last_transcript is not None) and (last_transcript.start < current_transcript.start):

            # Process VCF lines
            while (current_pos < current_transcript.start and line != ""):
                # check the allele frequency
                AF_pass = min_af <= 0
                if ';AF=' in line:
                    AF = float(line.split('AF=')[1].split(';')[0])
                    AF_pass = AF >= min_af

                # check validity of alleles
                val_pass = True
                REF, ALT = line.split(maxsplit=5)[3:5]
                if ((re.match(r'[CGTA]*[^CGTA]+[CGTA]*', REF)) or (re.match(r'[CGTA]*[^CGTA]+[CGTA]*', ALT))):
                    val_pass = False

                # check all transcripts in the queue
                if AF_pass and val_pass:
                    for transcript_entry in transcript_queue:

                        # check if the snp belongs to any of the exons
                        for exon in transcript_entry['exons']:
                            if (exon.start < current_pos):
                                if (exon.end > current_pos):
                                    transcript_entry['file_content'] += line
                                    break
                            else:
                                break   # exon starts after the mutation -> continue to another transcript

                vcf_linecount += 1
                line = vcf_file.readline()
                if line == "":
                    break

                current_pos = int(line.split(maxsplit=2)[1])
                vcf_id = line.split(maxsplit=3)[2]

                if (vcf_id == '.'):
                    # add an identifier = line cound
                    line = '\t'.join(line.split(maxsplit=2)[:2]) + '\t' + hex(vcf_linecount)[2:] + '\t' + line.split(maxsplit=3)[3]

            # remove passed transcripts from queue
            while (len(transcript_queue) > 0 and transcript_queue[0]['end'] < current_pos):
                df = pd.read_csv(StringIO(VCF_header + transcript_queue[0]['file_content']), sep='\t')
                df.to_csv(tmp_dir + '/' + transcript_queue[0]['ID'] + '.tsv', sep='\t', index=False, header=True)
                #result_dfs[transcript_queue[0]['ID']] = df
                transcript_queue.pop(0)

        # add the new transcript to the queue    
        exons = [ exon for exon in annotations_db.children(current_transcript, featuretype='exon', order_by='start') ]
        queue_entry = { 'transcript_obj': current_transcript, 'ID': current_transcript.id, 'exons': exons, 'start': current_transcript.start, 'end': current_transcript.end, 'file_content': "" }
        nearest_idx = bisect.bisect_left(KeyWrapper(transcript_queue, key=lambda x: x['end']), queue_entry['end'])
        transcript_queue.insert(nearest_idx, queue_entry)

        last_transcript = current_transcript

    # write the output for the remaining transcripts
    while (len(transcript_queue) > 0):
        df = pd.read_csv(StringIO(VCF_header + transcript_queue[0]['file_content']), sep='\t')
        df.to_csv(tmp_dir + '/' + transcript_queue[0]['ID'] + '.tsv', sep='\t', index=False, header=True)
        #result_dfs[transcript_queue[0]['ID']] = df
        transcript_queue.pop(0)

    return list(df.columns.values)