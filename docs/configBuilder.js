
function setGeneralParams(ensemblRelease, custom_tr_list, only_MANE, transcript_biotypes, contam_fasta_path, final_fasta_path, simplify_headers) {
    const general_config_html = "<p>\
        # ---------------- Global parameters ----------------<br>" +
        "Ensembl_FTP_URL: ftp.ensembl.org/pub/release-" + ensemblRelease + "/<br>" +
        "annotationFilename: Homo_sapiens.GRCh38." + ensemblRelease + ".chr_patch_hapl_scaff<br>" +
        "ensembl_release: " + ensemblRelease + "<br>" +
        "<br>" +
        "custom_transcript_list: \"" + custom_tr_list + "\"<br>"+
        "included_transcript_biotypes: \"" + transcript_biotypes + "\"<br>" +
        "only_MANE_select: " + (only_MANE ? "True" :  "False") + "<br>" +
        "<br>" +
        "contaminants_fasta: " + (contam_fasta_path.length === 0 ? "MISSING" : "\"" + contam_fasta_path + "\"") + "<br>" +
        "<br>" +
        "final_fasta_file: " + (final_fasta_path.length === 0 ? "MISSING" :  "\"" + final_fasta_path + "\"") + "<br>" +
        "fasta_simplify_headers: " + (simplify_headers ? "True" :  "False") + "<br>" +
        "<\p>"

    return general_config_html
}

function setProVarParams(use_ProVar, VCF_data, var_require_start, var_fasta_file, var_table_file, add_existing_haplo, haplo_added_fasta) {
    let vcf_datasets = ""

    VCF_data.forEach(element => {
        vcf_datasets += '&nbsp;&nbsp;' + element[0].replace(' ', '_') + ': { file: \"' + element[1] + "\", fasta_accession_prefix: " + element[0].replace(' ', '_') + '_' + ", min_af: " + element[2] + " }<br>"
    });

    const provar_config_html = "<p>\
        # ---------------- ProVar parameters ----------------<br>" +
        "use_ProVar: " + (use_ProVar ? "True" : "False") + '<br>' +
        (VCF_data.length === 0 ? "variant_vcf: {}<br>" : "variant_vcf: <br>" + vcf_datasets) +
        "<br>" + 
        "working_dir_name_var: variants_tmp<br><br>" +
        "var_require_start: " + (var_require_start ? "1" : "0") + '<br>' +
        "<br>" +
        "var_fasta_file: " + ((var_fasta_file.length == 0 && use_ProVar) ? "MISSING" :  "\"" + var_fasta_file + "\"") + "<br>" +
        "var_table_file: " + ((var_table_file.length == 0 && use_ProVar) ? "MISSING" :  "\"" + var_table_file + "\"") + "<br>" +
        "<br>" +
        "add_existing_haplo: " + (add_existing_haplo ? "True" : "False") + '<br>' +
        "haplo_added_fasta: " + ((haplo_added_fasta.length == 0 && add_existing_haplo) ? "MISSING" :  "\"" + haplo_added_fasta + "\"") + "<br>" +
        "</p>"

    return provar_config_html
}

function setProHapParams(use_ProHap, dataset_url, dataset_path, vcf_filename, samples_filename, freq_threshold, count_threshold, MAF_threshold, MAF_field, x_par1_to, x_par2_from, haplo_require_start, haplo_ignore_UTR, haplo_skip_start_lost, haplo_fasta_file, haplo_table_file) {
    if ((dataset_path.length > 0) && !dataset_path.endsWith('/')) {
        dataset_path = dataset_path + '/'
    }
    
    const prohap_config_html = "<p>\
    # ---------------- ProHap parameters ----------------<br>" +
    "use_ProHap: " + (use_ProHap ? "True" : "False") + '<br>' +
    "<br>" +
    "phased_FTP_URL: " + (dataset_url.length == 0 ? "\"\"" :  "\"" + dataset_url + "\"") + '<br>' +
    "phased_local_path: " + (dataset_path.length == 0 ? "\"\"" :  "\"" + dataset_path + "\"") + '<br>' +
    "phased_vcf_file_name: " + (vcf_filename.length == 0 ? "MISSING" :  "\"" + vcf_filename + "\"") + '<br>' +
    "sample_metadata_file: " + (samples_filename.length == 0 ? "MISSING" :  "\"" + samples_filename + "\"") + '<br>' +
    "<br>" +
    "working_dir_name_haplo: haplotypes_tmp<br><br>" +
    "phased_min_af: " + MAF_threshold + "<br>" +
    "phased_af_field: " + MAF_field + "<br>" +
    "haplo_min_freq: " + freq_threshold + "<br>" +
    "haplo_min_count: " + count_threshold + "<br>" +
    "<br>" +    
    "x_par1_to: " + x_par1_to + "<br>" +
    "x_par2_from: " + x_par2_from + "<br>" +
    "<br>" +
    "haplo_require_start: " + (haplo_require_start ? "1" : "0") + "<br>" +
    "haplo_ignore_UTR: " + (haplo_ignore_UTR ? "1" : "0") + "<br>" +
    "haplo_skip_start_lost: " + (haplo_skip_start_lost ? "1" : "0") + "<br>" +
    "max_cores: 3<br>" +
    "<br>" +
    "haplo_fasta_file: " + ((haplo_fasta_file.length == 0 && use_ProHap) ? "MISSING" :  "\"" + haplo_fasta_file + "\"") + "<br>" +
    "haplo_table_file: " + ((haplo_table_file.length == 0 && use_ProHap) ? "MISSING" :  "\"" + haplo_table_file + "\"") + "<br>" +
    "</p>"

    return prohap_config_html
}