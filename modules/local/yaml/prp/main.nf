process create_prp_yaml {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), val(lims_id), val(sample_name), path(quast), path(postalignqc), path(mlst), path(chewbbaca), path(amrfinder), path(resistance), path(resfinder_meta), path(serotypefinder), path(serotypefinder_meta), path(virulence), path(virulencefinder_meta), path(shigapass), path(emmtyper), path(bam), path(bai), path(nextflow_run_info), path(vcf), path(mykrobe), path(tbprofiler), path(bracken), path(bracken), path(bracken)
    path reference_genome
    path reference_genome_idx
    path reference_genome_gff

  output:
    tuple val(sample_id), path(output), emit: yaml

  script:
    output = "${sample_id}_prp.yaml"
    def access_dir = params.symlink_dir ? params.symlink_dir : params.outdir
    def amrfinder_kv = amrfinder ? "'amrfinder': '${amrfinder}'," : ""
    def bracken_kv = bracken ? "'kraken': '${bracken}'," : ""
    def bam_kv = bam ? "'uri': '${access_dir}/${params.species_dir}/${params.bam_dir}/${bam}'," : ""
    def bai_kv = bai ? "'index_uri': '${access_dir}/${params.species_dir}/${params.bam_dir}/${bai}'," : ""
    def chewbbaca_kv = chewbbaca ? "'chewbbaca': '${chewbbaca}'," : ""
    def emmtyper_kv = emmtyper ? "'emmtyper': '${emmtyper}'," : ""
    def lims_id_kv = lims_id ? "'lims_id': '${lims_id}'," : ""
    def mlst_kv = mlst ? "'mlst': '${mlst}'," : ""
    def mykrobe_kv = mykrobe ? "'mykrobe': '${mykrobe}'," : ""
    def postalignqc_kv = postalignqc ? "'postalnqc': '${postalignqc}'," : "" 
    def quast_kv = quast ? "'quast': '${quast}'," : ""
    def reference_genome_kv = reference_genome ? "'ref_genome_sequence': '${reference_genome}'," : ""
    def reference_genome_gff_kv = reference_genome_gff ? "'ref_genome_annotation': '${reference_genome_gff}'," : ""
    def resfinder_kv = resistance ? "'resfinder': '${resistance}'," : ""
    def resfinder_meta_lv = resfinder_meta ? "'${resfinder_meta}, " : ""
    def nextflow_run_info_kv = nextflow_run_info ? "'nextflow_run_info': '${nextflow_run_info}'," : ""
    def serotypefinder_kv = serotypefinder ? "'serotypefinder': '${serotypefinder}'," : ""
    def serotypefinder_meta_lv = serotypefinder_meta ? "'${serotypefinder_meta}, " : ""
    def shigapass_kv = shigapass ? "'shigapass': '${shigapass}'," : ""
    def ska_dir = task.ext.args_ska ?: ""
    def ska_index_kv = ska_index ? "'ska_index': '${ska_dir}/${ska_index}'," : ""
    def sourmash_dir = task.ext.args_sourmash ?: ""
    def sourmash_signature_kv = sourmash_signature ? "'sourmash_signature': '${ska_dir}/${sourmash_signature}'," : ""
    def tb_grading_rules_bed_array = tb_grading_rules_bed ? "{'name': 'tbdb grading rules bed', 'type': 'bed', 'uri': '${tb_grading_rules_bed}'}," : ""
    def tbdb_bed_array = tbdb_bed ? "{'name': 'tbdb bed', 'type': 'bed', 'uri': '${tbdb_bed}'}," : ""
    def tbprofiler_kv = tbprofiler ? "'tbprofiler': '${tbprofiler}'," : ""
    def vcf_kv = vcf ? "'uri': '${access_dir}/${params.species_dir}/${params.vcfDir}/${vcf}'," : ""
    def virulencefinder_kv = virulence ? "'virulencefinder': '${virulence}'," : ""
    def virulencefinder_meta_lv = virulencefinder_meta ?: ""
    """
    #!/usr/bin/env python
    import yaml

    # Define the prp_input
    prp_input = {
        ${amrfinder_kv}
        ${chewbbaca_kv}
        "groups": [${params.species_dir}],
        "igv_annotations": [
            {
                "name": "Cool variants",
                "type": "variant",
                ${vcf_kv}
            },
            {
                "name": "Read coverage",
                "type": "alignment",
                ${bam_kv}
                ${bai_kv}
            },
            ${tb_grading_rules_bed_array}
            ${tbdb_bed_array}
        ],
        ${bracken_kv}
        ${emmtyper_kv}
        ${lims_id_kv}
        ${mlst_kv}
        ${mykrobe_kv}
        ${nextflow_run_info_kv}
        ${postalnqc_kv}
        ${quast_kv}
        ${reference_genome_kv}
        ${reference_genome_gff_kv}
        ${resfinder_kv}
        "sample_id": ${sample_id},
        "sample_name": ${sample_name},
        ${serotypefinder_kv}
        ${shigapass_kv}
        ${ska_index_kv}
        "software_info": [ ${resfinder_meta_lv}${serotypefinder_meta_lv}${virulencefinder_meta_lv} ],
        ${sourmash_signature_kv}
        ${tbprofiler_kv}
        ${virulencefinder_kv}
    }

    with open("${output}", "w") as fout:
        yaml.dump(prp_input, fout, default_flow_style=False)
    """

  stub:
    output = "${sample_id}_prp.yaml"
    """
    touch ${output}
    """
}
