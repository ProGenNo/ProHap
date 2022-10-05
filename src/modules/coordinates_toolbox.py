from Bio.Seq import Seq
from numpy import ceil, floor

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


# check if we have an alteration of the start codon (either inframe indel before it, or stop loss)
# return new start location, -1 if start lost
def check_start_change(original_start, original_rf, variant_rna_loc, ref_len, alt_len, ignore_frameshift):
    if (variant_rna_loc < original_start+3):
        if (variant_rna_loc + ref_len > original_start):
            return -1, -1   # original start codon affected by change

        if (abs(alt_len - ref_len) % 3) != 0: # frameshift before start codon
            if (ignore_frameshift):
                return original_start + (alt_len - ref_len), (original_rf  + (alt_len - ref_len)) % 3
            return -1, -1    

        # stop codon might be shifted by inframe indel
        return original_start + (alt_len - ref_len), original_rf

    # change happening after start codon
    return original_start, original_rf

def get_affected_codons(cdna, allele_loc, allele_len, reading_frame, protein_start):
    alleles_protein = []        # residues directly affected (ignoring prossible frameshift), stored in a list for all three reading frames
    protein_location = []       # location of these residues in the protein (can be negative if in 5' UTR), creating a list as it can differ with reading frame
    
    if (reading_frame == -1):
        for rf in range(3):                    
            protein_location.append(int(floor((allele_loc - rf) / 3)))
    else:
        protein_location = [int(floor((allele_loc - reading_frame) / 3) -  protein_start)]

    bpFrom = int(floor((allele_loc - max(reading_frame, 0)) / 3) * 3 + max(reading_frame, 0))   # if reading frame is unknown, assume 0 and add other reading frames later
    bpFrom = max(max(bpFrom, 0), reading_frame)                                                 # in case the beginning of the change is before the reading frame start

    bpTo = int(ceil((allele_loc + allele_len - max(reading_frame, 0)) / 3) * 3 + max(reading_frame, 0))

    if (bpTo-bpFrom > 2): # make sure we have at least 1 codon covered
        affected_codons = Seq(cdna[bpFrom:bpTo])
        alleles_protein = [str(affected_codons.transcribe().translate())]
    else:
        alleles_protein = ['-']

    if reading_frame == -1:
        for rf in [1,2]:
            bpFrom = int(floor((allele_loc - rf) / 3) * 3 + rf) 
            bpFrom = max(max(bpFrom, 0), rf)                                                    
            bpTo = int(ceil((allele_loc + allele_len - rf) / 3) * 3 + rf)

            if (bpTo-bpFrom > 2): # make sure we have at least 1 codon covered
                affected_codons = Seq(cdna[bpFrom:bpTo])
                alleles_protein.append(str(affected_codons.transcribe().translate()))
            else:
                alleles_protein.append('-')

    return alleles_protein, protein_location