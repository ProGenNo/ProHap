
'''
Computes the position of the mutation in the RNA sequence. 
Checks whether the reference allele intersects a splice junction - truncates the sequences (both reference and alternative, if applicable) if so.
Special case: Mutation reaches over an intron into another exon. Probably covered but not tested.
'''
def get_rna_position(transcript_id, dna_location, ref_allele, alt_allele, exons):
    ref_len = len(ref_allele)
    alt_len = len(alt_allele)
    rna_location = 0
    found = False
    mutation_intersects_intron = None

    # find the corresponding exon - see how many nucleotides were there before
    for exon_idx,exon in enumerate(exons):
        # exon before the mutation -> remember the full length
        if (exon.end < dna_location):
            rna_location += (exon.end - exon.start + 1)
        
        # exon where the mutation happens -> remember position within
        # check for allele sequences ovrlapping borders of the exon
        elif (exon.start < (dna_location + ref_len)):
            
            # check if the mutation covers the intron before
            if (exon.start > dna_location):
                intronic_len = exon.start - dna_location
                ref_allele = ref_allele[intronic_len:]
                alt_allele = alt_allele[intronic_len:]

                ref_len = ref_len - intronic_len
                alt_len = len(alt_allele)

                dna_location += intronic_len

                mutation_intersects_intron = exon_idx

            rna_location += (dna_location - exon.start)
            found = True

            if (dna_location + ref_len > exon.end):
                remaining_length = exon.end - dna_location + 1
                mutation_intersects_intron = exon_idx + 1

                # check if the mutation does not reach into the next exon
                if exon_idx < (len(exons) - 1) and (dna_location + ref_len > exons[exon_idx+1].start):
                    next_exon = exons[exon_idx+1]
                    start_again = next_exon.start - dna_location

                    ref_allele = ref_allele[:remaining_length] + ref_allele[start_again:]
                    ref_len = len(ref_allele)
                else:
                    ref_allele = ref_allele[:remaining_length]
                    ref_len = remaining_length

                # if there is an insertion that prolongs an exon, keep it, 
                # only truncate the alternative allele if the reference overlaps
                if (dna_location + alt_len > exon.end):
                    remaining_length = exon.end - dna_location + 1

                    # check if the mutation does not reach into the next exon
                    if exon_idx < (len(exons) - 1) and (dna_location + alt_len > exons[exon_idx+1].start):
                        next_exon = exons[exon_idx+1]
                        start_again = next_exon.start - dna_location

                        alt_allele = alt_allele[:remaining_length] + alt_allele[start_again:]
                        alt_len = len(alt_allele)
                    else:
                        alt_allele = alt_allele[:remaining_length]
                        alt_len = remaining_length

            # remember if we change the last or first 3 letters in the exon
            elif (exon.end - dna_location + ref_len < 3):
                mutation_intersects_intron = exon_idx + 1
            
            elif (dna_location - exon.start < 3):
                mutation_intersects_intron = exon_idx

            break

    if not found:
        print(transcript_id + ': DNA location ' + str(dna_location) + ' is not in an exon.')
    
    return rna_location, ref_allele, ref_len, alt_allele, alt_len, mutation_intersects_intron


def get_rna_position_simple(transcript_id,dna_location, exons):
    rna_location = 0
    found = False

    # find the corresponding exon - see how many nucleotides were there before
    for exon in exons:
        # exon before the location -> remember the full length
        if (exon.end < dna_location):
            rna_location += (exon.end - exon.start + 1)
        elif (exon.start <= dna_location):
            found = True
            rna_location += (dna_location - exon.start)
            break

    if not found:
        raise Exception(transcript_id + ': DNA location ' + str(dna_location) + ' is not in an exon.')

    return rna_location