process spatyper {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path(output), emit: tsv 
    tuple val("${task.process}"), val('spatyper'), eval("spaTyper --version 2>&1 | sed 's/spaTyper //'"), topic: versions, emit: versions_spatyper

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    output = "${sample_id}_spatyper.tsv"
    """
    spaTyper -f ${assembly} --output ${output} ${args}
    """

    stub:
    output = "${sample_id}_spatyper.tsv"
    """
    touch ${output}
    """
}
