process sourmash {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path(output), emit: signature
    tuple val("${task.process}"), val('sourmash'), eval("sourmash --version 2>&1 | sed 's/^.*sourmash // ; s/ .*//'"), topic: versions, emit: versions_sourmash

    script:
    def args = task.ext.args ?: ''
    output = "${sample_id}.sig"
    """
    sourmash sketch dna ${args} ${assembly} -o ${output}
    """

    stub:
    output = "${sample_id}.sig"
    """
    touch ${output}
    """
}
