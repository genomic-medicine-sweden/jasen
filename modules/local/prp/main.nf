process format_jasen {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(yaml)

    output:
    tuple val(sample_id), path(output), emit: json
    path "*versions.yml"              , emit: versions

    script:
    output = "${sample_id}_result.json"
    """
    prp format-jasen \\
      --sample ${yaml} \\
      --output ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     prp:
      version: \$(echo \$(prp --version 2>&1) | sed 's/prp, version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

    stub:
    output = "${sample_id}_result.json"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     prp:
      version: \$(echo \$(prp --version 2>&1) | sed 's/prp, version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}

process format_cdm {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(yaml)

    output:
    tuple val(sample_id), path(output), emit: json

    script:
    output = "${sample_id}_qc_result.json"
    """
    prp format-cdm \\
      --sample ${yaml} \\
      --output ${output}
    """

    stub:
    output = "${sample_id}_qc_result.json"
    """
    touch ${output}
    """
}
