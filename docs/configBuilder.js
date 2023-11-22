
function setGeneralParams(ensemblRelease, custom_tr_list, transcript_biotypes, final_fasta_path) {
    const general_config_html = "<p>\
        # ---------------- Global parameters ----------------<br>" +
        "Ensembl_FTP_URL: ftp.ensembl.org/pub/release-" + ensemblRelease + "/<br>" +
        "annotationFilename: Homo_sapiens.GRCh38." + ensemblRelease + ".chr_patch_hapl_scaff<br>" +
        "ensembl_release: " + ensemblRelease + "<br>" +
        "<br>" +
        "custom_transcript_list: \"" + custom_tr_list + "\"<br>"+
        "included_transcript_biotypes: \"" + transcript_biotypes + "\"<br>" +
        "<br>" +
        "final_fasta_file: " + (final_fasta_path.length == 0 ? "MISSING" :  "\"" + final_fasta_path + "\"") + "<br>" +
        "<\p>"

    return general_config_html
}

function setProVarParams(include_variants, VCF_data, var_require_start, var_fasta_file, var_table_file) {
    let vcf_datasets = ""

    VCF_data.forEach(element => {
        vcf_datasets += '&nbsp;&nbsp;' + element[0].replace(' ', '_') + ': { file: \"' + element[1] + "\", fasta_accession_prefix: " + element[0].replace(' ', '_') + '_' + ", min_af: " + element[2] + " }<br>"
    });

    const provar_config_html = "<p>\
        # ---------------- ProVar parameters ----------------<br>" +
        "include_var: " + (include_variants ? "True" : "False") + '<br>' +
        (VCF_data.length === 0 ? "variant_vcf: {}<br>" : "variant_vcf: <br>" + vcf_datasets) +
        "<br>" + 
        "working_dir_name_var: variants_tmp<br><br>" +
        "var_require_start: " + (var_require_start ? "1" : "0") + '<br>' +
        "<br>" +
        "var_fasta_file: " + ((var_fasta_file.length == 0 && include_variants) ? "MISSING" :  "\"" + var_fasta_file + "\"") + "<br>" +
        "var_table_file: " + ((var_table_file.length == 0 && include_variants) ? "MISSING" :  "\"" + var_table_file + "\"") + "<br>" +
        "</p>"

    return provar_config_html
}

function setProHapParams(include_haplotypes, dataset_url, dataset_path, vcf_filename, freq_threshold, count_threshold, MAF_threshold, MAF_field, haplo_require_start, haplo_ignore_UTR, haplo_skip_start_lost, haplo_fasta_file, haplo_table_file) {
    if ((dataset_path.length > 0) && !dataset_path.endsWith('/')) {
        dataset_path = dataset_path + '/'
    }
    
    const prohap_config_html = "<p>\
    # ---------------- ProHap parameters ----------------<br>" +
    "include_haplo: " + (include_haplotypes ? "True" : "False") + '<br>' +
    "<br>" +
    "phased_FTP_URL: " + (dataset_url.length == 0 ? "\"\"" :  "\"" + dataset_url + "\"") + '<br>' +
    "phased_local_path: " + (dataset_path.length == 0 ? "\"\"" :  "\"" + dataset_path + "\"") + '<br>' +
    "phased_vcf_file_name: " + (vcf_filename.length == 0 ? "MISSING" :  "\"" + vcf_filename + "\"") + '<br>' +
    "<br>" +
    "working_dir_name_haplo: haplotypes_tmp<br><br>" +
    "phased_min_af: " + MAF_threshold + "<br>" +
    "phased_af_field: " + MAF_field + "<br>" +
    "haplo_min_freq: " + freq_threshold + "<br>" +
    "haplo_min_count: " + count_threshold + "<br>" +
    "<br>" +
    "haplo_require_start: " + (haplo_require_start ? "1" : "0") + "<br>" +
    "haplo_ignore_UTR: " + (haplo_ignore_UTR ? "1" : "0") + "<br>" +
    "haplo_skip_start_lost: " + (haplo_skip_start_lost ? "1" : "0") + "<br>" +
    "max_cores: 3<br>" +
    "<br>" +
    "haplo_fasta_file: " + ((haplo_fasta_file.length == 0 && include_haplotypes) ? "MISSING" :  "\"" + haplo_fasta_file + "\"") + "<br>" +
    "haplo_table_file: " + ((haplo_table_file.length == 0 && include_haplotypes) ? "MISSING" :  "\"" + haplo_table_file + "\"") + "<br>" +
    "</p>"

    return prohap_config_html
}