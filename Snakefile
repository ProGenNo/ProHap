configfile: "config.yaml"

CHROMOSOMES = [str(x) for x in list(range(1, 23))] + ['X']

rule all:
    input:
        final_fasta="results/haplotypes_nc/haplo_all.fa"

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
        db="data/gtf/" + config['annotationFilename'] + "_chr{chr}.db",
        tr="data/chr{chr}_transcripts_noncoding.txt"
    shell:
        "python3 src/parse_gtf.py -i {input} -o {output.db} -transcript_list {output.tr}"

rule compute_haplotypes:
    input:
        db="data/gtf/" + config['annotationFilename'] + "_chr{chr}.db",
        tr="data/chr{chr}_transcripts_noncoding.txt",
        vcf="data/1000genomes_GRCh38_vcf/ALL.chr{chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf",
        fasta="data/fasta/total_cdnas.fa",
        samples="igsr_samples.tsv"
    output:
        csv="results/haplotypes_nc/haplo_chr{chr}.tsv",
        fasta="results/haplotypes_nc/haplo_chr{chr}.fa"
    params:
        log_file="log/chr{chr}.log"
    threads: 8
    shell:
        "python3 src/prohap.py "
        "-i {input.vcf} -db {input.db} -transcripts {input.tr} -cdna {input.fasta} -s {input.samples} "
        "-chr {wildcards.chr} -af 0.01 -foo 0.01 -acc_prefix enshap_{wildcards.chr} -id_prefix haplo_chr{wildcards.chr} "
        "-threads 8 -log {params.log_file} -output_csv {output.csv} -output_fasta {output.fasta} "

rule merge_fasta:
    input:
        expand("results/haplotypes_nc/haplo_chr{chr}.fa", chr=CHROMOSOMES)
    output:
        "results/haplotypes_nc/haplo_all.fa"
    params:
        input_file_list = ' '.join(expand("results/haplotypes_nc/haplo_chr{chr}.fa", chr=CHROMOSOMES))
    shell:
        "cat {params.input_file_list} > {output}"