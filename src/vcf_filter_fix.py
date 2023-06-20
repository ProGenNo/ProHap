import argparse

parser = argparse.ArgumentParser(
        description='Reads the VCF file, parses multi-allelic variants into multiple lines, and filters out variants under the MAF threshold.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input VCF")

parser.add_argument("-af", dest="min_af", required=False, type=float,
                    help="Allele Frequency (AF) lower threshold - default 0", default=0)

parser.add_argument("-o", dest="output_file", required=True,
                    help="output VCF")

args = parser.parse_args()

def get_MAF(info):
    if ';AF=' in info:
        return info.split(';AF=')[1].split(';')[0]
    elif ';MAF=' in info:
        return info.split(';MAF=')[1].split(';')[0]
    elif '\tAF=' in info:
        return info.split('\tAF=')[1].split(';')[0]
    elif '\tMAF=' in info:
        return info.split('\tMAF=')[1].split(';')[0]

    return "-1"


# read the header of the VCF - keep only the last line of the header
VCF_header = ""

vcf_file = open(args.input_file, 'r')
outfile = open(args.output_file, 'w')

line = vcf_file.readline()

while (line != "" and line.startswith('#')):
    VCF_header += line
    line = vcf_file.readline()

outfile.write(VCF_header)

# check if the VCF has any valid lines
if (line == ''):
    outfile.close()
    vcf_file.close()
    exit()

while (line != ""):
    ALT = line.split(maxsplit=5)[4]
    INFO = line.split(maxsplit=8)[7]
    MAF = get_MAF(INFO)

    if (',' in ALT):
        CHR = line.split(maxsplit=1)[0]
        POS = line.split(maxsplit=1)[0]
        ID = line.split(maxsplit=1)[0]
        REF = line.split(maxsplit=1)[0]

        for i,allele in enumerate(ALT.split(',')):
            allele_maf = float(MAF.split(',')[i])

            if (allele_maf >= args.min_af):
                invalid_gts = [1,2,3]
                invalid_gts.remove(i+1)
                GTs = line.split(maxsplit=9)[-1]
                for gt_id in invalid_gts:
                    GTs = GTs.replace(str(gt_id) + '|', '0|')
                    GTs = GTs.replace('|' + str(gt_id), '|0')

                GTs = GTs.replace(str(i+1) + '|', '1|')
                GTs = GTs.replace('|' + str(i+1), '|1')                

                outfile.write('\t'.join(line.split(maxsplit=4)[:-1] + [allele] + line.split(maxsplit=7)[5:-1]))
                outfile.write('\tMAF=' + str(allele_maf) + '\tGT\t' + GTs)

    else:
        allele_maf = float(MAF)
        if (allele_maf >= args.min_af):
            outfile.write(line)

    line = vcf_file.readline()

vcf_file.close()
outfile.close()