# ProHap & ProVar
Proteogenomics database-generation tool for protein haplotypes and variants 

## Input & Usage
Required ingredients:
 - GTF annotation file (Ensembl - downloaded automatically by Snakemake)
 - cDNA FASTA file (Ensembl - downloaded automatically by Snakemake)
 - (optional) ncRNA FASTA file (Ensembl - downloaded automatically by Snakemake)
 - ProHap: VCF with phased genotpyes, one file per chromosome \(such as [1000 Genomes Project](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/) - downloaded automatically by Snakemake\)
 - ProVar: VCF, single file per dataset

Required software: [Snakemake](https://snakemake.readthedocs.io/en/stable/) & [Conda](https://docs.conda.io/en/latest/)

Usage:
 1. Create a configuration file called `config.yaml` based on the instructions in `config_file_example`
 2. Test Snakemake with a dry-run: `snakemake -c<# provided cores> -n -q`
 2. Run the Snakemake pipeline to create your protein database: `snakemake -c<# provided cores> -p --use-conda`

### Example of the configuration file:
Below is an example of the `config.yaml` file to run ProHap on the 1000 Genomes Project on GRCh38 (1kGP) data set and two custom VCF files with the following parameters:
 - Minor allele frequency (MAF) threshold for 1kGP variants: 0.01
 - MAF threshold for custom VCF files: 0 (i.e. no threshold)
 - Haplotype frequency threshold: 0.005
 - Ensembl version: 108
 - Include only transcripts specified in `data/transcripts_reference_108.csv`
 - Ignore transcripts without mORF annotation for 1kGP and custom VCF files
 - Keep UTR variants when computing the list of haplotypes, but remove the UTR sequences from the final database
 - Use 3 CPU cores per variant
 
```
1kGP_FTP_URL: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/
1kGP_vcf_file_name: "ALL.chr{chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf"

Ensembl_FTP_URL: ftp.ensembl.org/pub/release-108/
annotationFilename: Homo_sapiens.GRCh38.108.chr_patch_hapl_scaff

working_dir_name_haplo: haplotypes_tmp
working_dir_name_var: variants_tmp

variant_vcf:
   RareVar: { file: "data/clinvar_vcf/selected_clinvar.vcf", fasta_accession_prefix: "clinvar_", min_af: 0 }
   InhouseVar: { file: "data/inhouse_vcf/inhouse_variation.csv", fasta_accession_prefix: "rarevar_", min_af: 0 }

custom_transcript_list: "data/transcripts_reference_108.csv"
included_transcript_biotypes: "all"

include_haplo: True
include_var: False

haplo_min_freq: 0.005
haplo_min_count: 0
1kGP_min_af: 0.01

var_require_start: 1
haplo_require_start: 1
haplo_ignore_UTR: 0
haplo_skip_start_lost: 1

max_cores: 3

final_fasta_file: "results/final_db.fa"
haplo_fasta_file: "results/haplotypes.fa"
var_fasta_file: "results/vcf_variants.fa"
haplo_table_file: "results/haplotypes.tsv"
var_table_file: "results/vcf_variants.tsv"
```

## Analyzing the peptide-spectrum matches
Once you obtain a list of peptide-spectrum matches (PSMs), you can use a script provided in this repository \([peptides_annotate_variation.py](https://github.com/ProGenNo/ProHap/blob/main/src/analysis/peptides_annotate_variation.py)\) to map the peptides back to the respective protein haplotype / variant sequences, and map the identified variants back to their genetic origin. An example of the required input file (tab-separated):

```
PSMId Sequence Proteins Positions
psm_1 ALPCGHCPEEWITYSNSCYYIGK prot_f5d7;prot_22054;prot_3712e;prot_3ac94;prot_4055f 110;110;110;110;110
psm_2 LGCVLMAWALYLSLGVLWVAQMLLELFPAPILR prot_b40;prot_272f;prot_27a25;prot_2c28b;prot_2dac8 30;2;19;2;30
psm_3 ISVGVAGDLNTVTMK prot_b40;prot_3c88;prot_9b74;prot_9cf6;prot_b882 15;4;15;4;4
```

The annotation script can be used with the following parameters:

```
python src/analysis/peptides_annotate_variation.py 
    -i <input file of PSMs> 
    -hap <table of haplotypes produced by ProHap, if haplotypes are included in the search database> 
    -var <table of variants produced by ProVar, if individual variants are included in the search database> 
    -f <optimized FASTA file provided by ProHap/ProVar> 
    -tr_id data/protein_transcript_ids_108.csv <protein and transcript ID mapping - provided in this repository for Ensembl v.108> 
    -g_id data/gene_transcript_ids_108.csv <transcript and gene ID mapping - provided in this repository for Ensembl v.108> 
    -t <# threads> 
    -o <output filename> 
    -log <log file>
```

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

### Metadata table
Metadata file provided in a tab-separated text-file format. The columns are:
 - `chromosome`
 - `TranscriptID`
 - `transcript_biotype`: Biotype of the matching transcript in Ensembl.
 - `HaplotypeID`: ID of the haplotype sequence, matching to the ID in the FASTA entry description.
 - `VCF_IDs`: IDs of the matching lines in the VCF file of provided
 - `DNA_changes`: List of changes in the format POS:REF>ALT, mapped to the DNA coordinates within the chromosome
 - `allele_frequencies`: List of allele vrequencies of the variants in cluded in this haplotype
 - `cDNA_changes`: List of changes in the format POS:REF>ALT, mapped to the coordinates within the cDNA of this transcript
 - `all_protein_changes`: List of changes in the format POS:REF>ALT, mapped to the coordinates within the protein sequence. The start codon is at position 0, so if a change happens in the 5' untranslated region (UTR), its coordinates within the protein are negative.
 - `protein_changes`: List of changes in the protein excluding synonymous mutations.
 - `reading_frame`: Canonical reading frame for this transcript, if known.
 - `protein_prefix_length`: Number of codons in the 5' UTR
 - `splice_sites_affected`: List of splice sites affected by a mutation, if any. (Splice site 0 happens between exon 1 and 2)
 - `occurrence_count`: Number of occurrences of this haplotype within the participants of the 1000 Genomes project (or within the individuals provided in the phased genotype VCF)
 - `frequency`: Frequency of this haplotype within the participants of the 1000 Genomes project (or within the individuals provided in the phased genotype VCF)
