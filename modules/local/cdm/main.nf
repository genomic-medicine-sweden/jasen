process export_to_cdm {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(qc), val(sequencing_run), val(lims_id), val(sample_name)
    val species

    output:
    tuple val(sample_id), path(output), emit: cdm

    script:
    output = "${sample_id}.cdmpy"
    sequencing_run = sequencing_run ? "--sequencing-run ${sequencing_run}" : ""
    lims_id = lims_id ? "--lims-id ${lims_id}" : ""
    def sample_name_arg = "'" + sample_name.replace("'", "'\"'\"'") + "'"
    """
    cat <<'END_CDM' > ${output}
${sequencing_run} --sample-id ${sample_name_arg} --assay ${species} --qc ${params.outdir}/${params.species_dir}/${params.cdm_dir}/${qc} ${lims_id}
END_CDM
    """

    stub:
    output = "${sample_id}.cdmpy"
    """
    touch ${output}
    """
}
