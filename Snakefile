configfile: "config.yaml"

CHROMOSOMES = [str(x) for x in list(range(1, 23))] + ['X']

rule all:
        input:
                in1=expand("data/gtf/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff_chr{chr}.gtf", chr=CHROMOSOMES),
                #in2=expand("data/chr{chr}/ready", chr=CHROMOSOMES)

rule download_vcf:
        output:
                "data/1000genomes_GRCh38_vcf/ALL.chr{chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf"
        shell:
                "wget " + config['1000GsURL'] + "ALL.chr{wildcards.chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -O {output}.gz  && gunzip {output}.gz"

rule download_gtf:
        output:
                "data/gtf/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf"
        shell:
                "wget ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf.gz -O {output}.gz && gunzip {output}.gz; "

# filter the GTF so that only features on one chromosome are present:
rule split_gtf:
        input:
            in1="data/gtf/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf"
        output:
            "data/gtf/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff_chr{chr}.gtf"
        shell:
            "grep \"^#\" {input} > {output}; "
            "grep \"^{wildcards.chr}\s\" {input} >> {output}"

# create the DB files from GTF for each chromosome
rule parse_gtf:
        input:
            "data/gtf/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff_chr{chr}.gtf"
        output:
            "data/gtf/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff_chr{chr}.db"
        shell:
            "python3 src/parse_gtf.py -i {input} -o {output}"

rule fragment_vcf:
        input:
                vcf = "data/1000genomes_GRCh38_vcf/ALL.chr{chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf",
                gene_list = "data/gene_list/chr{chr}.csv"
        output:
                out_dummy=temp("data/chr{chr}/ready")
        shell:
                "python3 src/fragment_vcf.py -i {input.vcf} -g {input.gene_list} -d data/chr{wildcards.chr}"

