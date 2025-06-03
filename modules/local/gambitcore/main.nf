process gambitcore {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(fasta)
    path gambit_db

    output:
    tuple val(sample_id), path(output), emit: tsv
    path "*versions.yml"              , emit: versions

    script:
    def args = task.ext.args ?: ''
    output = "${sample_id}_gambitcore.tsv"
    """
    gambitcore \\
        ${gambit_db} \\
        ${fasta} \\
        ${args} \\
        > ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     gambitcore:
      version: \$(echo \$(gambitcore --version 2>&1) | sed 's/gambit, version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

    stub:
    output = "${sample_id}_gambitcore.tsv"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     gambitcore:
      version: \$(echo \$(gambitcore --version 2>&1) | sed 's/gambit, version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
