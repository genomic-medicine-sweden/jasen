process bracken {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(report)
    path database

    output:
    tuple val(sample_id), path(output)       , emit: output
    tuple val(sample_id), path(output_report), emit: report
    tuple val("${task.process}"), val('bracken'), val('2.8'), topic: versions, emit: versions_bracken

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    output = "${sample_id}_bracken.out"
    output_report = "${sample_id}_bracken.report"
    """
    bracken \\
    ${args} \\
    -d ${database} \\
    -i ${report} \\
    -o ${output} \\
    -w ${output_report}
    """

    stub:
    output = "${sample_id}_bracken.out"
    output_report = "${sample_id}_bracken.report"
    """
    touch ${output}
    touch ${output_report}
    """
}
