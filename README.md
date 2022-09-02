# ProHap & ProVar
Proteogenomics database-generation tool for protein haplotypes and variants 

## Input & Usage
Required ingredients:
 - GTF annotation file (Ensembl)
 - cDNA FASTA file (Ensembl)
 - (optional) ncRNA FASTA file (Ensembl)
 - contaminant sequences FASTA (such as https://www.thegpm.org/crap/)
 - ProHap: VCF with phased genotpyes, one file per chromosome \(such as [1000 Genomes Project](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/)\)
 - ProVar: VCF, one file per chromosome

Required software: snakemake & conda

Usage:
 1. create a configuration file called `config.yaml` based on the instructions in `config_file_example`
 2. run the Snakemake pipeline to create your protein database: `snakemake -c<# provided cores> -p --use-conda`

## Output
### FASTA protein database
The protein sequences are first split into sub-sequences by start and stop codon positions, and then duplicate sequences are aggregated into one FASTA entry. The resulting file has the following format:
```
>tag|accession|<positions_within_protein> <protein_IDs> <protein_starts> <matching_proteins> <reading_frames>
PROTEINSEQUENCE
```
Possible tag values are:
 - `generic_cont`: At least one of the matching sequences belongs to a contaminant.
 - `generic_ref`: No matching contaminant, at least one of the matching sequences belongs to a canonical protein.
 - `generic_var`: No matching contaminant or canonical protein, at least one of the matching sequences belongs to a variant protein.
 - `generic_hap`: No matching contaminant, canonical or variant protein, all of the matching sequences belong to a non-canonical protein haplotype.

The tag values for haplotypes and variants are customizable in the config file. 

The fields included in the description of the FASTA elements are the following:
 - `positions_within_protein`: position of matching sub-sequences within the whole protein sequence, delimited by semicolon
 - `protein_IDs`: IDs of the sub-sequences after splitting the whole protein (redundant)
 - `protein_starts`: positions of the start residue (usually M) within the whole protein, if known (0 otherwise)
 - `matching_proteins`: IDs of the whole protein sequences matching to this sub-sequence. Variant and haplotype IDs can be mapped to the metadata table provided.
 - `reading_frames`: Reading frames in which the matching proteins are translated, if known.
