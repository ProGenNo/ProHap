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

Using ProHap with the full 1000 Genomes Project data set (as per default) requires about 1TB disk space!

Usage:
 1. Clone this repository: `git clone https://github.com/ProGenNo/ProHap.git; cd ProHap/;`
 2. Create a configuration file called `config.yaml` based on the instructions in `config_file_example`, or using https://progenno.github.io/ProHap/
 3. Test Snakemake with a dry-run: `snakemake -c<# provided cores> -n -q`
 4. Run the Snakemake pipeline to create your protein database: `snakemake -c<# provided cores> -p --use-conda`

### Example: ProHap on 1000 Genomes
In the first usage example, we provide a small example dataset taken from the 1000 Genomes Project on GRCh38. We will ProHap to create a database of protein haplotypes aligned with Ensembl v.111 (January 2024).

This example was tested with Ubuntu 22.04.3 LTS. Windows users are encouraged to use the [Windows Subsystem for Linux](https://ubuntu.com/desktop/wsl) to run the ProHap / ProVar pipeline. 

Requirements: Install Conda / Mamba and Snakemake using [this guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba). Hardware requirements: 4 CPU cores, ~5 GB disk space, < 5 GB RAM.
  
Use the following commands to run this example:

```
# Clone this repository:
git clone https://github.com/ProGenNo/ProHap.git ;
cd ProHap;

# Unpack the sample dataset
cd sample_data ;
gunzip sample_1kGP_common_global.tar.gz ;
tar xf sample_1kGP_common_global.tar ;
cd .. ;

# Copy the configuration to config.yaml
cp config_example1.yaml config.yaml ;

# Activate the snakemake conda environment and run the pipeline
conda activate snakemake ;
snakemake --cores 4 -p --use-conda ;
```

## Using the database for proteomic searches
Once you obtain a list of peptide-spectrum matches (PSMs), you can use a pipeline provided in this repository \([peptide_annotation](https://github.com/ProGenNo/ProHap/tree/main/peptide_annotation)\) to map the peptides back to the respective protein haplotype / variant sequences, and map the identified variants back to their genetic origin. For the usage and details, please refer to the following [wiki page](https://github.com/ProGenNo/ProHap/wiki/Using-the-database-for-proteomic-searches).

## Output
The ProHap / ProVar pipeline produces three kinds of output files. Below is a brief description, please refer to the [wiki page](https://github.com/ProGenNo/ProHap/wiki/Output-files) for further details.

1. *Concatenated FASTA file*: The main result of the pipeline is the concatenated FASTA file, consisting of the ProHap and/or ProVar output, reference sequences from Ensembl, and common contaminant sequences \([cRAP](https://www.thegpm.org/crap/)\). The file can be used with any search engine, but is optimized for compatibility with [SearchGUI](http://compomics.github.io/projects/searchgui) and [PeptideShaker](http://compomics.github.io/projects/peptide-shaker).
2. *Metadata table*: Additional information on the variant / haplotype sequences produced by the pipeline, such as genomic coordinates of the variants covered, variant consequence type, etc.
3. *cDNA translations FASTA*: FASTA file contains the original translations of variant / haplotype cDNA sequences prior to any optimization, the removal of UTR sequences, and merging with canonical proteins and contaminants.
