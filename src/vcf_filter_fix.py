import argparse
import gzip
import glob

parser = argparse.ArgumentParser(
        description='Reads the VCF file, parses multi-allelic variants into multiple lines, and filters out variants under the MAF threshold.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input VCF")

parser.add_argument("-af_field", dest="af_field", required=False, type=str,
                    help="Allele Frequency (AF) field name - default \"AF\"", default="AF")

parser.add_argument("-af", dest="min_af", required=False, type=float,
                    help="Allele Frequency (AF) lower threshold - default 0", default=0)

parser.add_argument("-chr", dest="chromosome", required=True,
                    help="chromosome being processed (e.g., 1, 12 or X)")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output VCF")

args = parser.parse_args()

def get_MAF(info):
    if ';' + args.af_field + '=' in info:
        return info.split(';' + args.af_field + '=')[1].split(';')[0]
    elif '\t' + args.af_field + '=' in info:
        return info.split('\t' + args.af_field + '=')[1].split(';')[0]
    elif info.startswith(args.af_field + '='):
        return info.split(args.af_field + '=')[1].split(';')[0]

    return "-1"

# if there are multiple files for the same chromosome (specified, e.g., as "chr22_*.vcf"), process all of them in the alphabetical order
input_files = sorted(glob.glob(args.input_file))

# open the output file, keep it gzipped
outfile = gzip.open(args.output_file, 'wt')

# Remember the number of original entries and the number of filtered ones
total_VCF_entries = 0
valid_VCF_entries = 0

for i in range(len(input_files)):

    # read the header of the VCF
    VCF_header = ""

    if input_files[i].endswith('.gz'):
        vcf_file = gzip.open(input_files[i], 'rt')
    else:
        vcf_file = open(input_files[i], 'r')

    line = vcf_file.readline()

    while (line != "" and line.startswith('#')):
        VCF_header += line
        line = vcf_file.readline()

    # only write the header once, assuming all the other VCFs will follow the same format
    if (i == 0):
        outfile.write(VCF_header)

    # check if the VCF has any valid lines
    if (line == ''):
        vcf_file.close()
        print(input_files[i] + ' is empty!')
        continue

    while (line != ""):
        ALT = line.split(maxsplit=5)[4]
        INFO = line.split(maxsplit=8)[7]
        GT_IDX = line.split(maxsplit=9)[8].split(':').index('GT')
        
        # If not thresholding by MAF, ignore any possible formatting errors
        MAF = "0"
        if (args.min_af > 0):
            MAF = get_MAF(INFO)

        total_VCF_entries += 1

        # ProHap does not work with multi-allelic variants -> make a separate line in the output VCF for each allele
        if (',' in ALT):
            for i,allele in enumerate(ALT.split(',')):
                if (args.min_af > 0):
                    allele_maf = float(MAF.split(',')[i])
                else:
                    allele_maf = 0

                # keep only the i-th allele as the alternatvie (1), mark all the other alleles in the genotypes as 0
                if (allele_maf >= args.min_af):
                    invalid_gts = list(range(1,100))    # assuming that a variant will not have more than 100 possible alleles
                    invalid_gts.remove(i+1)

                    # format the string containing all the genotypes
                    GTs = [ GT.split(':')[GT_IDX] for GT in line.split()[9:] ]

                    # on sex chromosomes for male individuals, there might be just a single allele and the genotype won't be formated as "x|x" -> keep the valid allele as the first, and add a 0 to the genotype for consistency
                    GTs_str = '\t'.join([ ((GT + '|0') if (('|' not in GT) and ('/' not in GT)) else GT) for GT in GTs ])

                    # mark all the other alleles as 0
                    for gt_id in invalid_gts:
                        GTs_str = GTs_str.replace(str(gt_id) + '|', '0|')
                        GTs_str = GTs_str.replace('|' + str(gt_id), '|0')

                    # mark this allele as 1
                    GTs_str = GTs_str.replace(str(i+1) + '|', '1|')
                    GTs_str = GTs_str.replace('|' + str(i+1), '|1')                

                    outfile.write('\t'.join(line.split(maxsplit=4)[:-1] + [allele] + line.split(maxsplit=7)[5:-1]))
                    outfile.write('\tMAF=' + str(allele_maf) + '\tGT\t' + GTs_str + '\n')
                    valid_VCF_entries += 1

        else:
            allele_maf = float(MAF)
            if (allele_maf >= args.min_af):
                # format the string containing all the genotypes
                GTs = [ GT.split(':')[GT_IDX] for GT in line.split()[9:] ]

                # on sex chromosomes for male individuals, there might be just a single allele and the genotype won't be formated as "x|x" 
                # -> keep the valid allele as the first, and add a 0 to the genotype for consistency
                GTs_str = '\t'.join([ ((GT + '|0') if (('|' not in GT) and ('/' not in GT)) else GT) for GT in GTs ])

                outfile.write('\t'.join(line.split(maxsplit=8)[:-1]) + '\tGT\t' + GTs_str + '\n')
                valid_VCF_entries +=1

        line = vcf_file.readline()

    vcf_file.close()

outfile.close()

print(('Chr ' + args.chromosome + ':'), 'Original VCF(s):', total_VCF_entries, 'lines')
print(('Chr ' + args.chromosome + ':'), 'Filtered VCF:', valid_VCF_entries, 'lines')
