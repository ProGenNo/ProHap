configfile: "config.yaml"

VARIANT_VCF_FILES = config['variant_vcf']
CHROMOSOMES = config['chromosomes']

WORKING_DIR_NAME_HAPLO = config['working_dir_name_haplo']
WORKING_DIR_NAME_VAR = config['working_dir_name_var']

rule all:
    input:
        final_fasta=expand('{proxy}', proxy=[config['final_fasta_file']] if not config['fasta_simplify_headers'] else ['.'.join(config['final_fasta_file'].split('.')[:-1]) + '_simplified.fasta']),
        var_table=expand('{proxy}', proxy=[config['var_table_file']] if config["use_ProVar"] else []),
        haplo_table=expand('{proxy}', proxy=[config['haplo_table_file']] if config["use_ProHap"] else []),
        var_fasta=expand('{proxy}', proxy=[config['var_fasta_file']] if config["use_ProVar"] else []),
        haplo_fasta=expand('{proxy}', proxy=[config['haplo_fasta_file']] if config["use_ProHap"] else []),

rule download_vcf:
    output:
        temp("data/vcf/phased/" + config['phased_vcf_file_name'])
    shell:
        "mkdir -p data/vcf/phased ; "
        "wget " + config['phased_FTP_URL'] + config['phased_vcf_file_name'].replace('{chr}', '{wildcards.chr}') + " -O {output}"

rule download_gtf:
    output:
        temp("data/gtf/" + config['annotationFilename'] + ".gtf")
    shell:
        "mkdir -p data/gtf ; "
        "wget " + config['Ensembl_FTP_URL'] + "gtf/homo_sapiens/" + config['annotationFilename'] + ".gtf.gz -O {output}.gz && gunzip {output}.gz; "

rule parse_gtf_whole:
    input:
        "data/gtf/" + config['annotationFilename'] + ".gtf"
    output:
        temp("data/gtf/" + config['annotationFilename'] + ".db")
    conda: "envs/prohap.yaml"
    shell:
        "python3 src/parse_gtf.py -i {input} -o {output}"

rule get_transcript_list:
    input:
        "data/gtf/" + config['annotationFilename'] + ".db"
    output:
        temp("data/included_transcripts.csv")
    params:
        biotypes=config['included_transcript_biotypes']
    shell:
        "python3 src/get_transcript_list.py -i {input} -bio {params.biotypes} -o {output}"

rule download_cdnas_fasta:
    output:
        out1=temp("data/fasta/Homo_sapiens.GRCh38.ncrna.fa"),
        out2=temp("data/fasta/Homo_sapiens.GRCh38.cdna.all.fa")
    shell:
        "mkdir -p data/fasta ; "
        "wget " + config['Ensembl_FTP_URL'] + "fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz -O {output.out1}.gz && gunzip {output.out1}.gz; "
        "wget " + config['Ensembl_FTP_URL'] + "fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -O {output.out2}.gz && gunzip {output.out2}.gz; "

rule merge_cdnas_fasta:
	input:
		in1="data/fasta/Homo_sapiens.GRCh38.ncrna.fa",
		in2="data/fasta/Homo_sapiens.GRCh38.cdna.all.fa"
	output:
		"data/fasta/total_cdnas_" + str(config['ensembl_release']) + ".fa"
	shell:
		"cat {input.in1} > {output}; cat {input.in2} >> {output}"

rule download_reference_proteome:
    output:
        temp("data/fasta/Homo_sapiens.GRCh38.pep.all.fa")
    shell:
        "mkdir -p data/fasta ; "
        "wget " + config['Ensembl_FTP_URL'] + "fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz -O {output}.gz && gunzip {output}.gz; "

rule reference_filter_format:
    input:
        "data/fasta/Homo_sapiens.GRCh38.pep.all.fa"
    output:
        "data/fasta/ensembl_reference_proteinDB_" + str(config['ensembl_release']) + "_tagged.fa"
    conda: "envs/prohap.yaml"
    shell:
        "python3 src/fasta_format_headers.py -i {input} -o {output} -t _ensref -use_ENST 1 "

rule default_transcript_list:
    input:
        ref_fasta="data/fasta/ensembl_reference_proteinDB_" + str(config['ensembl_release']) + "_tagged.fa",
        annot="data/gtf/" + config['annotationFilename'] + ".db"
    output:
        temp("data/transcripts_reference_" + str(config['ensembl_release'])  + ".csv")
    params:
        MANE=int(config['only_MANE_select'])
    conda: "envs/prohap.yaml"
    shell:
        "python3 src/get_reference_ENST.py -i {input.ref_fasta} -annot {input.annot} -MANE {params.MANE} -o {output}"

rule reference_remove_stop:
    input:
        fasta="data/fasta/ensembl_reference_proteinDB_" + str(config['ensembl_release']) + "_tagged.fa",
        tr=expand('{proxy}', proxy=["data/transcripts_reference_" + str(config['ensembl_release']) + '.csv'] if config['only_MANE_select'] else [])
    output:
        temp("data/fasta/ensembl_reference_proteinDB_" + str(config['ensembl_release']) + "_clean.fa.gz")
    params:
        tr_filter=('-tr data/transcripts_reference_' + str(config['ensembl_release']) + '.csv') if config['only_MANE_select'] else ''
    conda: "envs/prohap.yaml"
    shell:
        "python3 src/remove_stop_codons.py -i {input.fasta} -o {output} -min_len 6 {params.tr_filter} "

rule contaminants_fix_headers:
    input:
        config['contaminants_fasta']
    output:
        "data/fasta/crap_tagged.fa.gz"
    conda: "envs/prohap.yaml"
    shell:
        "python3 src/fasta_format_headers.py -i {input} -o {output} -t _cont"

# filter the GTF so that only features on the desired chromosome are present
rule split_gtf:
    input:
        "data/gtf/" + config['annotationFilename'] + ".gtf"
    output:
        temp("data/gtf/" + config['annotationFilename'] + "_chr{chr}.gtf")
    shell:
        "grep \"^#\" {input} > {output}; "
        "grep \"^{wildcards.chr}\s\" {input} >> {output}"

# create the DB files from GTF for each chromosome
rule parse_gtf_chromosome:
    input:
        "data/gtf/" + config['annotationFilename'] + "_chr{chr}.gtf"
    output:
        temp("data/gtf/" + config['annotationFilename'] + "_chr{chr}.db")
    conda: "envs/prohap.yaml"
    shell:
        "python3 src/parse_gtf.py -i {input} -o {output}"

# ------------------------------------ ProVar rules ------------------------------------

rule split_variant_vcf:
    input:
        lambda wildcards: VARIANT_VCF_FILES[f"{wildcards.vcf}"]['file']
    output:
        temp("tmp/variants_{vcf}/ready")
    params:
        output_prefix="tmp/variants_{vcf}/variants"
    conda: "envs/prohap.yaml"
    shell:
        "python src/fragment_variant_vcf.py -i {input} -o {params.output_prefix} ; touch {output}"

rule compute_variants:
    input:
        db="data/gtf/" + config['annotationFilename'] + "_chr{chr}.db",
        tr=expand('{proxy}', proxy=[config['custom_transcript_list']] if len(config["custom_transcript_list"]) > 0 else ["data/included_transcripts.csv"]),
        fasta="data/fasta/total_cdnas_" + str(config['ensembl_release']) + ".fa",
        flag="tmp/variants_{vcf}/ready",
    output:
        tsv=temp("results/" + WORKING_DIR_NAME_VAR + "/variants_{vcf}/variants_chr{chr}.tsv.gz"),
        fasta=temp("results/" + WORKING_DIR_NAME_VAR + "/variants_{vcf}/variants_chr{chr}.fa" + ('.gz' if config['var_fasta_file'].endswith('.gz') else ""))
    params:
        input_vcf="tmp/variants_{vcf}/variants_chr{chr}.vcf.gz",
        output_cdna_file=("results/" + WORKING_DIR_NAME_VAR + "variants_cdna_chr{chr}.fa" + ('.gz' if config['var_cdna_file'].endswith('.gz') else "")) if (len(config['var_cdna_file']) > 0) else "",
        acc_prefix=lambda wildcards: VARIANT_VCF_FILES[f"{wildcards.vcf}"]['fasta_accession_prefix'],
        min_af=lambda wildcards: VARIANT_VCF_FILES[f"{wildcards.vcf}"]['min_af'],
        log_file="log/{vcf}_chr{chr}.log",
        #log_file="log/provar.log",
        tmp_dir="tmp/transcripts_{vcf}",
        require_start=config['var_require_start']
    conda: "envs/prohap.yaml"
    shell:
        "mkdir -p {params.tmp_dir}; mkdir -p log; mkdir -p results; "
        "python3 src/provar.py "
        "-i {params.input_vcf} -db {input.db} -transcripts {input.tr} -cdna {input.fasta} "
        "-chr {wildcards.chr} -acc_prefix {params.acc_prefix} -af {params.min_af} -require_start {params.require_start} "
        "-log {params.log_file} -tmp_dir {params.tmp_dir} -output_csv {output.tsv} -output_fasta {output.fasta} " + 
        ("-output_cdna_fasta \"{params.output_cdna_file}\" " if (len(config['var_cdna_file']) > 0) else "")

rule merge_var_tables_vcf:
    input:
        expand("results/" + WORKING_DIR_NAME_VAR + "/variants_{{vcf}}/variants_chr{chr}.tsv.gz", chr=CHROMOSOMES)
    output:
        temp("results/" + WORKING_DIR_NAME_VAR + "/variants_{vcf}/variants_all.tsv.gz")
    params:
        input_file_list = ','.join(expand("results/" + WORKING_DIR_NAME_VAR + "/variants_{{vcf}}/variants_chr{chr}.tsv.gz", chr=CHROMOSOMES))
    conda: "envs/prohap.yaml"
    shell:
        "python3 src/merge_tables.py -i {params.input_file_list} -o {output}"

rule merge_var_fasta_vcf:
    input:
        expand("results/" + WORKING_DIR_NAME_VAR + "/variants_{{vcf}}/variants_chr{chr}.fa" + ('.gz' if config['var_fasta_file'].endswith('.gz') else ""), chr=CHROMOSOMES)
    output:
        temp("results/" + WORKING_DIR_NAME_VAR + "/variants_{vcf}/variants_all.fa" + ('.gz' if config['var_fasta_file'].endswith('.gz') else ""))
    params:   
        input_file_list_cdna = expand("results/" + WORKING_DIR_NAME_VAR + "/variants_cdna_chr{chr}.fa" + ('.gz' if config['var_cdna_file'].endswith('.gz') else ""), chr=CHROMOSOMES),
        output_cdna_file = "results/" + WORKING_DIR_NAME_VAR + "/variants_{vcf}/variants_cdna_all.fa" + ('.gz' if config['var_cdna_file'].endswith('.gz') else "")
    shell:
        "cat {input} > {output}" + 
        ("cat {params.input_file_list_cdna} > {params.output_cdna_file} ; rm {params.input_file_list_cdna}" if (len(config['var_cdna_file']) > 0) else "")

rule merge_var_tables:
    input:
        expand("results/" + WORKING_DIR_NAME_VAR + "/variants_{vcf}/variants_all.tsv.gz", vcf=VARIANT_VCF_FILES.keys())
    output:
        config['var_table_file']
    params:
        input_file_list = ','.join(expand("results/" + WORKING_DIR_NAME_VAR + "/variants_{vcf}/variants_all.tsv.gz", vcf=VARIANT_VCF_FILES.keys()))
    conda: "envs/prohap.yaml"
    shell:
        "python3 src/merge_tables.py -i {params.input_file_list} -o {output}"

rule merge_var_fasta:
    input:
        expand("results/" + WORKING_DIR_NAME_VAR + "/variants_{vcf}/variants_all.fa" + ('.gz' if config['var_fasta_file'].endswith('.gz') else ""), vcf=VARIANT_VCF_FILES.keys())
    output:
        config['var_fasta_file']
    params:
        input_file_list_cdna = ','.join(expand("results/" + WORKING_DIR_NAME_VAR + "/variants_{vcf}/variants_cdna_all.fa" + ('.gz' if config['var_cdna_file'].endswith('.gz') else ""), vcf=VARIANT_VCF_FILES.keys())),
        output_cdna_file = config['var_cdna_file']
    shell:
        "cat {input} > {output}" + 
        ("cat {params.input_file_list_cdna} > {params.output_cdna_file} ; rm {params.input_file_list_cdna}" if (len(config['var_cdna_file']) > 0) else "")

rule var_fasta_remove_stop:
    input:
        config['var_fasta_file']
    output:
        temp("results/variants_all_clean.fa.gz")
    conda: "envs/prohap.yaml"
    shell:
        "python3 src/remove_stop_codons.py -i {input} -o {output} -min_len 6 "

# ------------------------------------ ProHap rules ------------------------------------

rule filter_phased_vcf:
    # Skip the input as it will be inferred by the script from the glob pattern
    #input:
    #    vcf=expand('{proxy}', proxy=[config['phased_local_path'] + config['phased_vcf_file_name']] if len(config["phased_local_path"]) > 0 else ["data/vcf/phased/" + config['phased_vcf_file_name']])
    output:
        temp("data/vcf/phased/chr{chr}_phased_filtered.vcf.gz")
    params:
        input_file_glob=(config['phased_local_path'] + config['phased_vcf_file_name']) if len(config["phased_local_path"]) > 0 else ("data/vcf/phased/" + config['phased_vcf_file_name']),
        AF_threshold=config['phased_min_af'],
        AF_field=config['phased_af_field']
    shell:
        "mkdir -p data/vcf/phased ; "
        "python3 src/vcf_filter_fix.py -i \'{params.input_file_glob}\' -chr {wildcards.chr} -af {params.AF_threshold} -af_field {params.AF_field} -o {output} "

rule compute_haplotypes:
    input:
        db="data/gtf/" + config['annotationFilename'] + "_chr{chr}.db",
        tr=expand('{proxy}', proxy=[config['custom_transcript_list']] if len(config["custom_transcript_list"]) > 0 else ["data/included_transcripts.csv"]),
        vcf="data/vcf/phased/chr{chr}_phased_filtered.vcf.gz",
        fasta="data/fasta/total_cdnas_" + str(config['ensembl_release']) + ".fa",
        samples=config['sample_metadata_file']
    output:
        csv=temp("results/" + WORKING_DIR_NAME_HAPLO + "/haplo_chr{chr}.tsv.gz"),
        fasta=temp("results/" + WORKING_DIR_NAME_HAPLO + "/haplo_chr{chr}.fa" + ('.gz' if config['haplo_fasta_file'].endswith('.gz') else "")),
    params:
        log_file="log/prohap_chr{chr}.log",
        tmp_dir="tmp/transcript_vcf_haplo",
        output_cdna_file=("results/" + WORKING_DIR_NAME_HAPLO + "/haplo_cdna_chr{chr}.fa" + ('.gz' if config['haplo_cdna_file'].endswith('.gz') else "")) if (len(config['haplo_cdna_file']) > 0) else "",
        require_start=config['haplo_require_start'],
        ignore_UTR=config['haplo_ignore_UTR'],
        skip_start_lost=config['haplo_skip_start_lost'],
        freq_threshold=config['haplo_min_freq'],
        count_threshold=config['haplo_min_count'],
        x_par1_to=config['x_par1_to'],
        x_par2_from=config['x_par2_from'],
        max_cores=config['max_cores']
    threads: config['max_cores']
    conda: "envs/prohap.yaml"
    shell:
        "mkdir -p {params.tmp_dir}; mkdir -p log; mkdir -p results; "
        "python3 src/prohap.py "
        "-i \"{input.vcf}\" -db \"{input.db}\" -transcripts \"{input.tr}\" -cdna \"{input.fasta}\" -s \"{input.samples}\" "
        "-chr {wildcards.chr} -min_hap_freq {params.freq_threshold} -min_hap_count {params.count_threshold} "
        "-acc_prefix enshap_{wildcards.chr} -id_prefix haplo_chr{wildcards.chr} -require_start {params.require_start} -ignore_UTR {params.ignore_UTR} -skip_start_lost {params.skip_start_lost} "
        "-x_par1_to {params.x_par1_to} -x_par2_from {params.x_par2_from} -threads {params.max_cores} -log \"{params.log_file}\" -tmp_dir \"{params.tmp_dir}\" -output_csv \"{output.csv}\" -output_fasta \"{output.fasta}\" " + 
        ("-output_cdna_fasta \"{params.output_cdna_file}\" " if (len(config['haplo_cdna_file']) > 0) else "")

rule merge_haplo_tables:
    input:
        expand("results/" + WORKING_DIR_NAME_HAPLO + "/haplo_chr{chr}.tsv.gz", chr=CHROMOSOMES)
    output:
        #config['haplo_table_file']
        temp("results/" + WORKING_DIR_NAME_HAPLO + "/haplo_all.tsv.gz")
    params:
        input_file_list = ','.join(expand("results/" + WORKING_DIR_NAME_HAPLO + "/haplo_chr{chr}.tsv.gz", chr=CHROMOSOMES))
    conda: "envs/prohap.yaml"
    shell:
        "python3 src/merge_tables.py -i {params.input_file_list} -o {output}"

rule extract_sample_names:
    input:
        "results/" + WORKING_DIR_NAME_HAPLO + "/haplo_all.tsv.gz"
    output:
        haplo_tsv=config['haplo_table_file'],
        samples=('.'.join(config['haplo_table_file'].split('.')[:-2]) + '_sampleIDs.tsv.gz') if config['haplo_table_file'].endswith('.gz') else ('.'.join(config['haplo_table_file'].split('.')[:-1]) + '_sampleIDs.tsv')
    conda: "envs/prohap.yaml"
    shell:
        "python3 src/haplo_extract_sample_names.py -hap_tsv {input} -o {output.haplo_tsv} -samples {output.samples} "    

rule merge_fasta:
    input:
        expand(("results/" + WORKING_DIR_NAME_HAPLO + "/haplo_chr{chr}.fa" + ('.gz' if config['haplo_fasta_file'].endswith('.gz') else "")), chr=CHROMOSOMES)
    output:
        config['haplo_fasta_file']
    params:
        input_file_list_cdna = expand("results/" + WORKING_DIR_NAME_HAPLO + "/haplo_cdna_chr{chr}.fa" + ('.gz' if config['haplo_cdna_file'].endswith('.gz') else ""), chr=CHROMOSOMES),
        output_cdna_file = config['haplo_cdna_file']
    shell:
        "cat {input} > {output} ; " + 
        ("cat {params.input_file_list_cdna} > {params.output_cdna_file} ; rm {params.input_file_list_cdna}" if (len(config['haplo_cdna_file']) > 0) else "")

rule haplo_fasta_remove_stop:
    input:
        config['haplo_fasta_file']
    output:
        temp("results/haplo_all_clean.fa.gz")
    conda: "envs/prohap.yaml"
    shell:
        "python3 src/remove_stop_codons.py -i {input} -o {output} -min_len 6 "

# ------------------------------------ post-processing rules ------------------------------------

rule added_fasta_remove_stop:
    input:
        config['haplo_added_fasta']
    output:
        temp("results/haplo_added_clean.fa.gz")
    conda: "envs/prohap.yaml"
    shell:
        "python3 src/remove_stop_codons.py -i {input} -o {output} -min_len 6 "

rule mix_with_reference_proteome:
    input:
        in1="data/fasta/ensembl_reference_proteinDB_" + str(config['ensembl_release']) + "_clean.fa.gz",
        in2="data/fasta/crap_tagged.fa.gz",
        in3=expand('{proxy}', proxy=["results/variants_all_clean.fa.gz"] if config["use_ProVar"] else []),
        in4=expand('{proxy}', proxy=["results/haplo_all_clean.fa.gz"] if config["use_ProHap"] else []),
        in5=expand('{proxy}', proxy=["results/haplo_added_clean.fa.gz"] if config["add_existing_haplo"] else []),
    output:
        temp("results/ref_contam_vcf_haplo_all_clean.fa.gz")		
    run:
        shell("cat {input.in1} {input.in2} > {output}; ")
        if config["use_ProVar"]:
            shell("cat {input.in3} >> {output}")
        if config["use_ProHap"]:
            shell("cat {input.in4} >> {output}")
        if config["add_existing_haplo"]:
            shell("cat {input.in5} >> {output}")

rule merge_duplicate_seq:
    input:
        "results/ref_contam_vcf_haplo_all_clean.fa.gz"
    output:
        temp("results/ref_contam_vcf_haplo_all_nodupl.fa.gz")
        #config['final_fasta_file']                         
    conda: "envs/prohap.yaml"
    shell:
        "python3 src/merge_duplicate_seq.py -i {input} -o {output} "

# UTRs in ProHap are removed by default (can be changed), but not in ProVar -> make sure all UTRs are removed
rule remove_UTR_seq:
    input:
        "results/ref_contam_vcf_haplo_all_nodupl.fa.gz"
    output:
        config['final_fasta_file']
    conda: "envs/prohap.yaml"
    shell:
        "python src/remove_UTR_seq.py -i {input} -o {output}"

rule simpify_fasta_headers:
    input:
        fasta=config['final_fasta_file'],
        # var_table=expand('{proxy}', proxy=[config['var_table_file']] if config["use_ProVar"] else []), ENST identifiers included in the protein IDs for ProVar - no need to search for them in the TSV
        haplo_table=expand('{proxy}', proxy=[config['haplo_table_file']] if config["use_ProHap"] else ([config['haplo_added_table']] if (config['use_ProVar'] and config['add_existing_haplo']) else [])),
        annot="data/gtf/" + config['annotationFilename'] + ".db"
    output:
        fasta=('.'.join(config['final_fasta_file'].split('.')[:-2]) + '_simplified.fasta.gz') if config['final_fasta_file'].endswith('.gz') else ('.'.join(config['final_fasta_file'].split('.')[:-1]) + '_simplified.fasta'),
        header=('.'.join(config['final_fasta_file'].split('.')[:-2]) + '_header.tsv.gz') if config['final_fasta_file'].endswith('.gz') else ('.'.join(config['final_fasta_file'].split('.')[:-1]) + '_header.tsv')
    conda: "envs/prohap.yaml"
    shell:
        "python src/fasta_simplify_headers.py -i {input.fasta} " +
        ("-hap_tsv {input.haplo_table} " if (config["use_ProHap"] or (config['use_ProVar'] and config['add_existing_haplo'])) else "") +
        "-db {input.annot} -o {output.fasta} -header {output.header}"