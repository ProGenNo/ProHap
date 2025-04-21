import pandas as pd
import os
import gzip

def check_open_file(parser, arg):
        if not os.path.exists(arg):
                parser.error("The file %s does not exist!" % arg)
        elif arg.endswith('.gz'):
                return gzip.open(arg, 'rt')  # return an open file handle
        else:
               return open(arg, 'r')

# returns an object containing all the sequences + metadata in the fasta file,  
# accessed by the stable element ID or accession in case of artificial identifier
# all_elements[elementID] = {'tag': tag, 'accession': accession, 'description': description, 'sequence': sequence}
def read_fasta(filename, truncate_accession = False):

    fasta_file = gzip.open(filename, 'rt') if (filename.endswith('.gz')) else open(filename, 'r')

    metadata = fasta_file.readline()    # line starting with '>'
    sequence = ""                       # all the other lines are considered sequences (considering also a multi-line format)

    all_elements = {}

    # read the reference sequence database
    while metadata != "":

        line = fasta_file.readline() 
        while not line.startswith('>') and not line == "":
            sequence += line[:-1]
            line = fasta_file.readline()

        tag = ''
        accession = ''
        description = ''

        if "|" in metadata:					# the header is at least partially formated
            metadata_parsed = metadata[1:].split('|')
            if 'generic' in metadata_parsed[0]:
                tag = metadata_parsed[0]
            else:
                tag = 'generic_' + metadata_parsed[0]

            if len(metadata_parsed) == 2:					# accession and description are potentially merged -> separate them
                if " " in metadata_parsed[1]:
                    accession = metadata_parsed[1].split(' ')[0]
                    description = metadata_parsed[1].split(' ', 1)[1]
                else:
                    accession = metadata_parsed[1]				# no description -> keep accession as it is
            elif len(metadata_parsed) == 3:					# descripton and accesson already separated -> keep
                accession = metadata_parsed[1]
                description = metadata_parsed[2]

        else:						                    # the header is not formated
            accession = metadata[1:].split(" ")[0]
            if " " in metadata:
                description = metadata.split(" ", 1)[1]

        if (truncate_accession):
            elementID = accession.split('.')[0]
        else:
            elementID = accession

        all_elements[elementID] = {'tag': tag, 'accession': accession, 'description': description.replace('\n', ''), 'sequence': sequence}

        metadata = line
        sequence = ""

    fasta_file.close()

    return all_elements

def check_vcf_df(in_df):
        result_columns = in_df.columns.values
        result_data = []
        for index,row in in_df.iterrows():
                ALT_alleles = row['ALT'].split(',')
                REF = row['REF']
                for ALT in ALT_alleles:
                        result_row = row.tolist()
                        result_row[4] = ALT
                        result_data.append(result_row)
        result_df = pd.DataFrame(columns=result_columns, data=result_data)
        return result_df

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
