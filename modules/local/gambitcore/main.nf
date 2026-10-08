process gambitcore {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(fasta)
    path gambit_db

    output:
    tuple val(sample_id), path(output), emit: tsv
    tuple val("${task.process}"), val('gambitcore'), eval("gambitcore --version 2>&1 | sed 's/^gambitcore //'"), topic: versions, emit: versions_gambitcore

    script:
    def args = task.ext.args ?: ''
    output = "${sample_id}_gambitcore.tsv"
    """
    gambitcore \\
        ${gambit_db} \\
        ${fasta} \\
        ${args} \\
        > ${output}
    """

    stub:
    output = "${sample_id}_gambitcore.tsv"
    """
    touch ${output}
    """
}
