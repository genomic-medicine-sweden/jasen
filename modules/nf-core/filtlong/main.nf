process filtlong {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path(output), emit: reads
    tuple val("${task.process}"), val('filtlong'), eval("filtlong --version 2>&1 | sed 's/Filtlong v//'"), topic: versions, emit: versions_filtlong

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    output = "${sample_id}_filtered.fastq.gz"
    """
    filtlong \\
        ${args} \\
        ${reads} \\
        2>| >(tee ${sample_id}_filtlong.log >&2) \\
        | gzip -n > ${output}
    """

    stub:
    output = "${sample_id}_filtered.fastq.gz"
    """
    touch ${output}
    """
}
