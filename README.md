# ProHap & ProVar
Proteogenomics database-generation tool for protein haplotypes and variants. Paper describing the tool: [![DOI](https://zenodo.org/badge/DOI/10.1038/s41592-024-02506-0.svg)](https://doi.org/10.1038/s41592-024-02506-0).

When using ProHap or ProVar in your publication, please cite: Vašíček, J., Kuznetsova, K.G., Skiadopoulou, D. et al. ProHap enables human proteomic database generation accounting for population diversity. _Nature Methods_ (2024). [https://doi.org/10.1038/s41592-024-02506-0](https://doi.org/10.1038/s41592-024-02506-0)

## Databases generated using ProHap
- Databases obtained from the common haplotypes of the 1000 Genomes Project along with metadata set can be found at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10149277.svg)](https://doi.org/10.5281/zenodo.10149277).
- Databases obtained from the common haplotypes of the Release 1.1 of the Haplotype Reference Consortium (HRC) along with metadata set can be found at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12671301.svg)](https://doi.org/10.5281/zenodo.12671301).
- Databases obtained from the preliminary release of the Human Pangenome Reference Consortium (HPRC) along with metadata set can be found at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12686818.svg)](https://doi.org/10.5281/zenodo.12686818).

> Note: The databases contain only common haplotypes (maf > 1 %), no individual-level data is available from the databases. For individual-level sequences, please run ProHap on the individual-level data.

## Input & Usage
Below is a brief overview, for details on input file format and configuration, please refer to the [Wiki page](https://github.com/ProGenNo/ProHap/wiki/Input-&-Usage).

Required input:
 - For ProHap: VCF with phased genotypes, one file per chromosome \(such as [1000 Genomes Project](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/) - downloaded automatically by Snakemake if URL is provided\)
 - For ProVar: VCF, single file per dataset. Multiple VCF files can be processed by ProVar in the same run.
 - FASTA file of contaminant sequences. These will then be added to the final FASTA, and tagged as contaminants. The default contaminant database is created by the [cRAP](https://www.thegpm.org/crap/) project, provided in this repository.
 - GTF annotation file (Ensembl - downloaded automatically by Snakemake)
 - cDNA FASTA file (Ensembl - downloaded automatically by Snakemake)
 - (optional) ncRNA FASTA file (Ensembl - downloaded automatically by Snakemake)

Required software: [Snakemake](https://snakemake.readthedocs.io/en/stable/) & [Conda](https://docs.conda.io/en/latest/). ProHap was tested with Ubuntu 22.04.3 LTS. Windows users are encouraged to use the [Windows Subsystem for Linux](https://ubuntu.com/desktop/wsl).

Using ProHap with the full 1000 Genomes Project data set (as per default) requires about 1TB disk space!

Usage:
 1. Clone this repository: `git clone https://github.com/ProGenNo/ProHap.git; cd ProHap/;`
 2. Create a configuration file called `config.yaml` using https://progenno.github.io/ProHap/. Please refer to the [Wiki page](https://github.com/ProGenNo/ProHap/wiki/Input-&-Usage) for details.
 3. Test Snakemake with a dry-run: `snakemake --cores <# provided cores> -n -q`
 4. Run the Snakemake pipeline to create your protein database: `snakemake --cores <# provided cores> -p --use-conda`

### Example: ProHap on 1000 Genomes
In the first usage example, we provide a small example dataset taken from the 1000 Genomes Project on GRCh38. We will use ProHap to create a database of protein haplotypes aligned with Ensembl v.115 (May 2025) using only MANE Select transcripts.

Expected runtime using 4 CPU cores: ~1 hour. Expected runtime using 23 CPU cores: ~30 minutes.

Requirements: Install Conda / Mamba and Snakemake using [this guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba). Minimum hardware requirements: 1 CPU core, ~5 GB disk space, 3 GB RAM.
  
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
Once you obtain a list of peptide-spectrum matches (PSMs), you can use a pipeline provided in the [PeptideAnnotator](https://github.com/ProGenNo/ProHap_PeptideAnnotator) repository to map the peptides back to the respective protein haplotype / variant sequences, and map the identified variants back to their genetic origin. For the usage and details, please refer to the following [wiki page](https://github.com/ProGenNo/ProHap/wiki/Using-the-database-for-proteomic-searches).

## Output
The ProHap / ProVar pipeline produces four kinds of output files. Below is a brief description, please refer to the [wiki page](https://github.com/ProGenNo/ProHap/wiki/Output-files) for further details.

1. *Concatenated FASTA file*: The main result of the pipeline is the concatenated FASTA file, consisting of the ProHap and/or ProVar output, reference sequences from Ensembl, and provided contaminant sequences. The file can be used with any search engine.
    * Optionally, headers are extracted and provided in an attached tab-separated file, and a gene name is added to each protein entry. 
2. *Metadata table*: Additional information on the variant / haplotype sequences produced by the pipeline, such as genomic coordinates of the variants covered, variant consequence type, etc.
3. *cDNA translations FASTA*: FASTA file contains the original translations of variant / haplotype cDNA sequences prior to any optimization, and merging with canonical proteins and contaminants. Note that if ignoring variation in UTRs (default configuration of ProHap), UTR sequences are not included here. Otherwise, they are kept in this file. 
4. *cDNA sequences FASTA*: Optional: FASTA file contains the variant / haplotype cDNA sequences prior to any optimization, and merging with canonical proteins and contaminants. Note that if ignoring variation in UTRs (default configuration of ProHap), UTR sequences are not included here. Otherwise, they are kept in this file. 

## Bug report and contribution
We welcome bug reports, suggestions of improvements, and contributions. Please do not hesitate to [open an issue](https://github.com/ProGenNo/ProHap/issues) or a [pull request](https://github.com/ProGenNo/ProHap/pulls).

## Code of Conduct
As part of our efforts toward delivering open and inclusive science, we follow the [Contributor Convenant](https://www.contributor-covenant.org/) [Code of Conduct for Open Source Projects](docs/CODE_OF_CONDUCT.md).

## Citation
When using ProHap and databases generated using ProHap, please cite the accompanying scientific publication [![DOI](https://zenodo.org/badge/DOI/10.1101/2023.12.24.572591.svg)](https://doi.org/10.1101/2023.12.24.572591).


