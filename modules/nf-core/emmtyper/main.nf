process emmtyper {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path(output), emit: tsv
    tuple val("${task.process}"), val('emmtyper'), eval("echo \$(emmtyper --version 2>&1) | sed -r 's/^.*emmtyper // ; s/ .*//'"), topic: versions, emit: versions_emmtyper

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    output = "${sample_id}_emmtyper.tsv"
    """
    emmtyper ${args} --output ${output} ${assembly}
    """

    stub:
    output = "${sample_id}_emmtyper.tsv"
    """
    touch ${output}
    """
}
