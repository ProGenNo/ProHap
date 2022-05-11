import gffutils
import argparse
import os.path
import bisect

# hepler class for bisect to sort objects
class KeyWrapper:
    def __init__(self, iterable, key):
        self.it = iterable
        self.key = key

    def __getitem__(self, i):
        return self.key(self.it[i])

    def __len__(self):
        return len(self.it)

    def insert(self, index, item):
        self.it.insert(index, item)


parser = argparse.ArgumentParser(
        description='Read a VCF file (sorted by position!) and creates a VCF file for each of the transcripts. Works only per chromosome.')

def is_valid_file(parser, arg):
        if not os.path.exists(arg):
                parser.error("The file %s does not exist!" % arg)
        else:
                return open(arg, 'r')  # return an open file handle

parser.add_argument("-i", dest="input_file", required=True,
                        help="input VCF file", metavar="FILE",
                        type=lambda x: is_valid_file(parser, x))

parser.add_argument("-db", dest="annotation_db", required=True,
                    help="DB file created by gffutils from GTF")

parser.add_argument("-d", dest="output_dir", required=True,
                    help="output directory")

parser.add_argument("-ch", dest="chrom", required=True,
                    help="chromosome number")

parser.add_argument("-foo", dest="min_foo", required=False, type=float,
                    help="Frequency of Occurance (FoO) lower threshold - default 0", default=0)

args = parser.parse_args()

# Load the annotations database
annotations_db = gffutils.FeatureDB(args.annotation_db)

# read the header of the VCF
VCF_header = ""

line = args.input_file.readline()
vcf_linecount = 0

while (line != "" and line.startswith('#')):
    VCF_header += line
    line = args.input_file.readline()
    vcf_linecount += 1

# browse the chromosome in a sweep-line approach - assumes that the VCF file is sorted!
# keep a list of transcripts that intersect the current position of the sweep line -> assign the VCF line to all of these transcripts

transcript_queue = []               # queue of transcript objects inc. the exons, sorted by end position, each element aggregates the VCF file contents
file_handles = {}                   # file handles accessed by transcript ID
current_pos = int(line.split()[1])  # position of the current VCF entry 
transcript_count = 0
last_transcript = None

# iterate through all the transcripts - add the first one to the queue (and all others starting at the same location), and of each other, check if there is a gap that can be filled in by VCF entries
# i.e., process all the VCF entries that lay before the current transcript -> update the queue
for current_transcript in annotations_db.features_of_type('transcript', order_by='start'):

    if (last_transcript is not None) and (last_transcript.start < current_transcript.start):

        # TODO: process VCF lines
        while (current_pos < current_transcript.start and line != ""):
            # check all transcripts in the queue
            for transcript_entry in transcript_queue:

                # check if the snp belongs to any of the exons
                for exon in transcript_entry['exons']:
                    if (exon.start < current_pos):
                        if (exon.end > current_pos):
                            transcript_entry['file_content'] += line
                            break
                    else:
                        break   # exon starts after the mutation -> continue to another transcript

            line = args.input_file.readline()
            if line == "":
                break
            current_pos = int(line.split()[1])

        # remove passed transcripts from queue
        while (len(transcript_queue) > 0 and transcript_queue[0]['end'] < current_pos):
            file = open(args.output_dir + '/' + transcript_queue[0]['ID'] + '.vcf', 'w')
            file.write(VCF_header)
            file.write(transcript_queue[0]['file_content'])
            file.close()
            transcript_queue.pop(0)
        
    # add the new transcript to the queue    
    exons = [ exon for exon in annotations_db.children(current_transcript, featuretype='exon', order_by='start') ]
    queue_entry = { 'transcript_obj': current_transcript, 'ID': current_transcript.id, 'exons': exons, 'start': current_transcript.start, 'end': current_transcript.end, 'file_content': "" }
    nearest_idx = bisect.bisect_left(KeyWrapper(transcript_queue, key=lambda x: x['end']), queue_entry['end'])
    transcript_queue.insert(nearest_idx, queue_entry)

    transcript_count += 1

file = open(args.output_dir + '/' + "ready", 'w')
file.write("Ready")
file.close()

args.input_file.close()
print('Done')