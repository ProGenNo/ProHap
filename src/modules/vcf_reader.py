import bisect
from modules.common import KeyWrapper
from io import StringIO
import pandas as pd
import re

def check_vcf_line_validity(line, min_af):
    # check the allele frequency
    AF_pass = min_af <= 0
    if ';AF=' in line:
        AF = float(line.split(';AF=')[1].split(';')[0])
        AF_pass = AF >= min_af
    elif ';MAF=' in line:
        AF = float(line.split(';MAF=')[1].split(';')[0])
        AF_pass = AF >= min_af

    # check validity of alleles
    val_pass = True
    REF, ALT = line.split(maxsplit=5)[3:5]
    if ((re.match(r'[CGTA]*[^CGTA]+[CGTA]*', REF) and REF != '-') or (re.match(r'[CGTA,]*[^CGTA,]+[CGTA,]*', ALT) and ALT != '-')):
        val_pass = False

    return AF_pass and val_pass

def add_variants_to_transcripts(vcf_file_line, vcf_file, vcf_linecount, transcript_queue, current_pos, current_transcript, VCF_header, min_af, tmp_dir, finalize):
    # Process VCF lines
    while ((current_pos < current_transcript.start or finalize) and vcf_file_line != ""):
        valid = check_vcf_line_validity(vcf_file_line, min_af)

        # check all transcripts in the queue
        if valid:
            for transcript_entry in transcript_queue:

                # check if the snp belongs to any of the exons
                for exon in transcript_entry['exons']:
                    if (exon.start <= current_pos):
                        if (exon.end >= current_pos):
                            transcript_entry['file_content'] += vcf_file_line
                            break
                    else:
                        break   # exon starts after the mutation -> continue to another transcript

        vcf_linecount += 1
        vcf_file_line = vcf_file.readline()
        if vcf_file_line == "":
            break

        current_pos = int(vcf_file_line.split(maxsplit=2)[1])
        vcf_id = vcf_file_line.split(maxsplit=3)[2]

        if (vcf_id == '.'):
            # add an identifier = line cound
            vcf_file_line = '\t'.join(vcf_file_line.split(maxsplit=2)[:2]) + '\t' + hex(vcf_linecount)[2:] + '\t' + vcf_file_line.split(maxsplit=3)[3]

    # remove passed transcripts from queue
    while (len(transcript_queue) > 0 and (transcript_queue[0]['end'] < current_pos or finalize)):
        df = pd.read_csv(StringIO(VCF_header + transcript_queue[0]['file_content']), sep='\t')
        df.to_csv(tmp_dir + '/' + transcript_queue[0]['ID'] + '.tsv', sep='\t', index=False, header=True)
        #result_dfs[transcript_queue[0]['ID']] = df
        transcript_queue.pop(0)

    return VCF_header[:-1].split('\t'), vcf_file_line, vcf_linecount, transcript_queue, current_pos

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

    # check if the VCF has any valid lines
    if (line == ''):
        return []

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

            colnames, line, vcf_linecount, transcript_queue, current_pos = add_variants_to_transcripts(line, vcf_file, vcf_linecount, transcript_queue, current_pos, current_transcript, VCF_header, min_af, tmp_dir, False)

        # add the new transcript to the queue    
        exons = [ exon for exon in annotations_db.children(current_transcript, featuretype='exon', order_by='start') ]
        queue_entry = { 'transcript_obj': current_transcript, 'ID': current_transcript.id, 'exons': exons, 'start': current_transcript.start, 'end': current_transcript.end, 'file_content': "" }
        nearest_idx = bisect.bisect_left(KeyWrapper(transcript_queue, key=lambda x: x['end']), queue_entry['end'])
        transcript_queue.insert(nearest_idx, queue_entry)

        last_transcript = current_transcript

    colnames, line, vcf_linecount, transcript_queue, current_pos = add_variants_to_transcripts(line, vcf_file, vcf_linecount, transcript_queue, current_pos, current_transcript, VCF_header, min_af, tmp_dir, True)

    return colnames
