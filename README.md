# ProHap & ProVar
Proteogenomics database-generation tool for protein haplotypes and variants. Preprint describing the tool: [doi.org/10.1101/2023.12.24.572591](https://doi.org/10.1101/2023.12.24.572591). 

A database created using ProHap on the 1000 Genomes Project data set can be found at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10149278.svg)](https://doi.org/10.5281/zenodo.10149278).

## Input & Usage
Below is a brief overview, for details on input file format and configuration, please refer to the [Wiki page](https://github.com/ProGenNo/ProHap/wiki/Input-&-Usage).

Required ingredients:
 - GTF annotation file (Ensembl - downloaded automatically by Snakemake)
 - cDNA FASTA file (Ensembl - downloaded automatically by Snakemake)
 - (optional) ncRNA FASTA file (Ensembl - downloaded automatically by Snakemake)
 - For ProHap: VCF with phased genotypes, one file per chromosome \(such as [1000 Genomes Project](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/) - downloaded automatically by Snakemake\)
 - For ProVar: VCF, single file per dataset

Required software: [Snakemake](https://snakemake.readthedocs.io/en/stable/) & [Conda](https://docs.conda.io/en/latest/)

Using ProHap with the 1000 Genomes Project data set (as per default) requires about 1TB disk space!

Usage:
 1. Clone this repository: `git clone https://github.com/ProGenNo/ProHap.git; cd ProHap/;`
 2. Create a configuration file called `config.yaml` based on the instructions in `config_file_example`, or using https://progenno.github.io/ProHap/
 3. Test Snakemake with a dry-run: `snakemake -c<# provided cores> -n -q`
 4. Run the Snakemake pipeline to create your protein database: `snakemake -c<# provided cores> -p --use-conda`

## Using the database for proteomic searches
Once you obtain a list of peptide-spectrum matches (PSMs), you can use a pipeline provided in this repository \([peptide_annotation](https://github.com/ProGenNo/ProHap/tree/main/peptide_annotation)\) to map the peptides back to the respective protein haplotype / variant sequences, and map the identified variants back to their genetic origin. For the usage and details, please refer to the following [wiki page](https://github.com/ProGenNo/ProHap/wiki/Using-the-database-for-proteomic-searches).

## Output
The ProHap / ProVar pipeline produces three kinds of output files. Below is a brief description, please refer to the [wiki page](https://github.com/ProGenNo/ProHap/wiki/Output-files) for further details.

1. *Concatenated FASTA file*: The main result of the pipeline is the concatenated FASTA file, consisting of the ProHap and/or ProVar output, reference sequences from Ensembl, and common contaminant sequences \([cRAP](https://www.thegpm.org/crap/)\). The file can be used with any search engine, but is optimized for compatibility with [SearchGUI](http://compomics.github.io/projects/searchgui) and [PeptideShaker](http://compomics.github.io/projects/peptide-shaker).
2. *Metadata table*: Additional information on the variant / haplotype sequences produced by the pipeline, such as genomic coordinates of the variants covered, variant consequence type, etc.
3. *cDNA translations FASTA*: FASTA file contains the original translations of variant / haplotype cDNA sequences prior to any optimization, the removal of UTR sequences, and merging with canonical proteins and contaminants.
