import argparse
import pandas as pd

parser = argparse.ArgumentParser(
    description='Merge all the haplotype tables into one TSV file.')

parser.add_argument("-i", dest="input_filenames", required=True,
                    help="input files, comma-separated list")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output TSV file", metavar="FILE",
                    type=lambda x: open(x, 'w'))

args = parser.parse_args()

inputs = args.input_filenames.split(",")
dataframes = []

for i,inputFilename in enumerate(inputs):
    chromosome = inputFilename.split('chr', 1)[1].split('.', 1)[0]
    df = pd.read_csv(inputFilename, sep='\t', header=0)
    df['chromosome'] = chromosome

result = pd.concat(dataframes)
result.to_csv(args.output_file, index=False, header=True, sep='\t')

args.output_file.close()
