configfile: "config.yaml"

CHROMOSOMES = [str(x) for x in list(range(1, 23))] + ['X']

rule all:
        input:
                in1=expand("data/gtf/" + config['annotationFilename'] + "_chr{chr}.gtf", chr=CHROMOSOMES),
                in2=expand("results/gene_haplotypes2/gene_haplo_chr{chr}.tsv", chr=CHROMOSOMES),
                in3="data/fasta/total_cdnas.fa"

rule download_vcf:
        output:
                "data/1000genomes_GRCh38_vcf/ALL.chr{chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf"
        shell:
                "wget " + config['1000GsURL'] + "ALL.chr{wildcards.chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -O {output}.gz  && gunzip {output}.gz"

rule download_gtf:
        output:
                "data/gtf/" + config['annotationFilename'] + ".gtf"
        shell:
                "wget " + config['EnsemblFTP'] + "gtf/homo_sapiens/" + config['annotationFilename'] + ".gtf.gz -O {output}.gz && gunzip {output}.gz; "

rule download_cdnas_fasta:
        output:
                out1="data/fasta/Homo_sapiens.GRCh38.ncrna.fa",
                out2="data/fasta/Homo_sapiens.GRCh38.cdna.all.fa"
        shell:
                "wget " + config['EnsemblFTP'] + "fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz -O {output.out1}.gz && gunzip {output.out1}.gz; "
                "wget " + config['EnsemblFTP'] + "fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -O {output.out2}.gz && gunzip {output.out2}.gz; "

rule merge_cdnas_fasta:
	input:
		in1="data/fasta/Homo_sapiens.GRCh38.ncrna.fa",
		in2="data/fasta/Homo_sapiens.GRCh38.cdna.all.fa"
	output:
		"data/fasta/total_cdnas.fa"
	shell:
		"cat {input.in1} > {output}; cat {input.in2} >> {output}"

# filter the GTF so that only features on one chromosome are present:
rule split_gtf:
        input:
            "data/gtf/" + config['annotationFilename'] + ".gtf"
        output:
            "data/gtf/" + config['annotationFilename'] + "_chr{chr}.gtf"
        shell:
            "grep \"^#\" {input} > {output}; "
            "grep \"^{wildcards.chr}\s\" {input} >> {output}"

# create the DB files from GTF for each chromosome
rule parse_gtf:
        input:
            "data/gtf/" + config['annotationFilename'] + "_chr{chr}.gtf"
        output:
            "data/gtf/" + config['annotationFilename'] + "_chr{chr}.db"
        shell:
            "python3 src/parse_gtf.py -i {input} -o {output}"

# make a separate VCF for each transcript - only coding variants with AF passing threshold
rule fragment_vcf:
        input:
                vcf="data/1000genomes_GRCh38_vcf/ALL.chr{chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf",
                db="data/gtf/" + config['annotationFilename'] + "_chr{chr}.db"
        output:
                out_dummy="data/chr{chr}/ready"
        shell:
                "python3 src/fragment_vcf.py -i {input.vcf} -db {input.db} -d data/chr{wildcards.chr} -foo 0.01"

# create the list of gene haplotypes for each chromosome
rule gene_haplotypes:
        input:
                in_dummy="data/chr{chr}/ready",
                db="data/gtf/" + config['annotationFilename'] + "_chr{chr}.db"
        output:
                "results/gene_haplotypes2/gene_haplo_chr{chr}.tsv"
        shell:
                "python3 src/get_haplotypes.py -d data/chr{wildcards.chr} -db {input.db} -o {output}"

# process the gene haplotypes into cDNA and protein haplotypes
rule protein_haplotypes:
        input:
                gene_csv="results/gene_haplotypes2/gene_haplo_chr{chr}.tsv",
                db="data/gtf/" + config['annotationFilename'] + "_chr{chr}.db",
                cdna="data/fasta/total_cdnas.fa"
        output:
                csv="results/protein_haplotypes2/haplo_chr{chr}.tsv",
                fasta="results/protein_haplotypes2/haplo_chr{chr}.fa
        shell:
                "python3 src/translate_haplotypes.py -i {input.gene_csv} -db {input.db} -cdna {input.cdna} -acc_prefix enshap_{wildcards.chr}_ -output_csv {output.csv} -output_fasta {output.fasta}"

