import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Reads a list of discoverable variants and a list of peptides, outputs a plot with number of variants and haplotypes per haplotype frequency.')

parser.add_argument("-v", dest="variant_list_1", required=True,
                    help="variant list (CSV)", metavar="FILENAME")

parser.add_argument("-v2", dest="variant_list_2", required=True,
                    help="variant list (CSV)", metavar="FILENAME")

parser.add_argument("-hap", dest="haplo_table_1", required=True,
                    help="haplotype table (TSV)", metavar="FILENAME")

parser.add_argument("-hap2", dest="haplo_table_2", required=True,
                    help="haplotype table (TSV)", metavar="FILENAME")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output PNG file", metavar="FILENAME")

args = parser.parse_args()

print ('Reading', args.variant_list_1)
variants_df_1 = pd.read_csv(args.variant_list_1)

print ('Reading', args.variant_list_2)
variants_df_2 = pd.read_csv(args.variant_list_2)

print ('Reading', args.haplo_table_1)
haplo_df_1 = pd.read_table(args.haplo_table_1)

print ('Reading', args.haplo_table_2)
haplo_df_2 = pd.read_table(args.haplo_table_2)

print ('Aggregating data')
x = np.linspace(0, 0.05, 150)
max_varcount = max(len(variants_df_1), len(variants_df_2))
max_hapcount = max(len(haplo_df_1), len(haplo_df_2))
y_var_1 = []
y_hap_1 = []

y_var_2 = []
y_hap_2 = []

for freq in x:
    y_var_1.append((len(variants_df_1[variants_df_1['max_haplotype_frequency'] >= freq]) / max_varcount) * 100)
    y_hap_1.append((len(haplo_df_1[haplo_df_1['frequency'] >= freq]) / max_hapcount) * 100)
    y_var_2.append((len(variants_df_2[variants_df_2['max_haplotype_frequency'] >= freq]) / max_varcount) * 100)
    y_hap_2.append((len(haplo_df_2[haplo_df_2['frequency'] >= freq]) / max_hapcount) * 100)

print ('Creating the plot')
plt.figure()
plt.plot(x, y_var_1, color='navy', label='variants')
plt.plot(x, y_var_2, color='navy', label='variants (UTRs ignored)', linestyle='--')
plt.plot(x, y_hap_1, color='red', label='haplotypes')
plt.plot(x, y_hap_2, color='red', label='haplotypes (UTRs ignored)', linestyle='--')

plt.xlabel("Haploptype Frequency Threshold")
plt.ylabel('% of included variants / haplotypes')
plt.xticks([0,0.005,0.01,0.015,0.02,0.03,0.04,0.05], rotation=30)
plt.grid()
plt.legend()
plt.tight_layout()

plt.savefig(args.output_file, dpi=300)
