import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Reads a list of discoverable variants and a list of peptides, outputs a plot with number of variants and haplotypes per haplotype frequency.')

parser.add_argument("-v", dest="variant_list", required=True,
                    help="variant list (CSV)", metavar="FILENAME")

parser.add_argument("-hap", dest="haplo_table", required=True,
                    help="haplotype table (TSV)", metavar="FILENAME")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output PNG file", metavar="FILENAME")

args = parser.parse_args()

print ('Reading', args.variant_list)
variants_df = pd.read_csv(args.variant_list)

print ('Reading', args.haplo_table)
haplo_df = pd.read_table(args.haplo_table)

print ('Aggregating data')
x = np.linspace(0, 0.05, 150)
y_var = []
y_hap = []

for freq in x:
    y_var.append((len(variants_df[variants_df['max_haplotype_frequency'] >= freq]) / len(variants_df)) * 100)
    y_hap.append((len(haplo_df[haplo_df['frequency'] >= freq]) / len(haplo_df)) * 100)

print ('Creating the plot')
plt.figure()
plt.plot(x, y_var, label='variants')
plt.plot(x, y_hap, color='red', label='haplotypes')

plt.xlabel("Haploptype Frequency Threshold")
plt.ylabel('% of included variants / haplotypes')
plt.grid()
plt.legend()

plt.savefig(args.output_file, dpi=300)
