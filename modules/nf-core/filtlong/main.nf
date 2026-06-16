process filtlong {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path(output), emit: reads
    path "*versions.yml",               emit: versions

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

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     filtlong:
      version: \$(filtlong --version | sed 's/Filtlong v//')
      container: ${task.container}
    END_VERSIONS
    """

    stub:
    output = "${sample_id}_filtered.fastq.gz"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     filtlong:
      version: \$(filtlong --version | sed 's/Filtlong v//')
      container: ${task.container}
    END_VERSIONS
    """
}
