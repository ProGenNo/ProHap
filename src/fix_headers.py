import argparse
import os.path

parser = argparse.ArgumentParser(description='Format the protein headers as follows: >generic[your tag]|[protein accession]|[protein description]. Creates a single-line fasta.')

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return open(arg, 'r')  # return an open file handle

parser.add_argument("-i", dest="input_file", required=True,
                    help="input FASTA file", metavar="FILE",
                    type=lambda x: is_valid_file(parser, x))

parser.add_argument("-o", dest="output_file", required=True,
                    help="output FASTA file", metavar="FILE",
                    type=lambda x: open(x, 'w'))

parser.add_argument("-t", dest="tag", required=False, default="",
		    help="custom tag for protein identification")

parser.add_argument("-use_ENST", dest="use_ENST", required=False, default=0,
		    help="replace the current protein accession with the ENST identifier")

args = parser.parse_args()

print('')
print("Reading file", args.input_file.name)
print("Formatting protein headers.")

metadata = args.input_file.readline()   # line starting with '>'
sequence = ""                           # all the other lines are considered protein sequences (considering also a multi-line format)
variant_count = 0		        # counter used to create unique identifiers for variants

while metadata != "":

    line = args.input_file.readline()
    while not line.startswith('>') and not line == "":
        sequence += line[:-1]
        line = args.input_file.readline()

    tag = ''
    accession = ''
    description = ''

    if "|" in metadata:					# the header is at least partially formated
        metadata_parsed = metadata[1:-1].split('|')
        if 'generic' in metadata_parsed[0]:
            tag = metadata_parsed[0]
        else:
            tag = 'generic_' + metadata_parsed[0]

        if len(metadata_parsed) == 2:					# accession and description are potentially merged -> separate them
            if " " in metadata_parsed[1]:
                description = metadata_parsed[1].split(' ', 1)[1]
                if (not args.use_ENST) or ('ENST' not in description):
                    accession = metadata_parsed[1].split(' ')[0]
                else:
                    accession = 'ENST' + description.split('ENST',1)[1].split(maxsplit=1)[0].split('.',1)[0]
            else:
                accession = metadata_parsed[1]				# no description -> keep accession as it is
        elif len(metadata_parsed) == 3:					# descripton and accesson already separated -> keep
            description = metadata_parsed[2]

            if (not args.use_ENST) or ('ENST' not in description):
                accession = metadata_parsed[1]
            else:
                accession = 'ENST' + description.split('ENST',1)[1].split(maxsplit=1)[0].split('.',1)[0]

    else:						# the header is not formated
        tag = 'generic' + args.tag					# add the keyword "generic" and a custom tag (empty if not provided)
        accession = metadata[1:-1].split()[0]
        
        if " " in metadata:
            description = metadata[:-1].split(" ", 1)[1]
            if (args.use_ENST) and ('ENST' in description):
                accession = 'ENST' + description.split('ENST',1)[1].split(maxsplit=1)[0].split('.',1)[0]

    if 'matching_proteins:' not in description:
        description = description + ' matching_proteins:' + accession + '\n'

    if len(description) > 0:
        new_header = '>' + '|'.join([tag, accession, description])
    else:
        new_header = '>' + '|'.join([tag, accession])

    args.output_file.write(new_header)

    if sequence.endswith('\n'): 
        args.output_file.write(sequence)    
    else:    
        args.output_file.write(sequence + '\n')    

    metadata = line
    sequence = ""

args.output_file.close()
args.input_file.close()
print('Done, output file:', args.output_file.name)
