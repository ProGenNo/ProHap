
# returns an object containing all the sequences + metadata in the fasta file,  
# accessed by the stable element ID or accession in case of artificial identifier
# all_elements[elementID] = {'tag': tag, 'accession': accession, 'description': description, 'sequence': sequence}
def read_fasta(filename):
    fasta_file = open(filename, 'r')

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
            if 'generic' in tag:
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

        elementID = accession.split('.')[0]

        all_elements[elementID] = {'tag': tag, 'accession': accession, 'description': description.replace('\n', ''), 'sequence': sequence}

        metadata = line
        sequence = ""

    fasta_file.close()

    return all_elements