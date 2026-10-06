process quast {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(assembly)
    path reference

    output:
    tuple val(sample_id), path(output), emit: tsv
    tuple val("${task.process}"), val('quast'), eval("echo \$(quast.py --version 2>&1) | sed 's/^.*QUAST v//'"), topic: versions, emit: versions_quast

    script:
    def args = task.ext.args ?: ''
    output = "${sample_id}_quast.tsv"
    reference_command = reference ? "-r ${reference}" : ""
    output_dir = "quast_outdir"
    """
    quast.py ${args} ${assembly} ${reference_command} -o ${output_dir} -t ${task.cpus}
    cp ${output_dir}/transposed_report.tsv ${output}
    """

    stub:
    output = "${sample_id}_quast.tsv"
    """
    touch ${output}
    """
}
