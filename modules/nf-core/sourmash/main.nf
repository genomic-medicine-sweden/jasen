process sourmash {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path(output), emit: signature
    path "*versions.yml"              , emit: versions

    script:
    def args = task.ext.args ?: ''
    output = "${sample_id}.sig"
    """
    sourmash sketch dna ${args} ${assembly} -o ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     sourmash:
      version: \$(echo \$(sourmash --version 2>&1) | sed 's/^.*sourmash // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

    stub:
    output = "${sample_id}.sig"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     sourmash:
      version: \$(echo \$(sourmash --version 2>&1) | sed 's/^.*sourmash // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
