process mykrobe {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path(output), emit: csv
    path "*versions.yml"              , emit: versions

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_reads_arg = reads.size() == 2 ? "${reads.join(' ')}" : "${reads[0]}"
    output = "${sample_id}_mykrobe.csv"
    """
    mykrobe predict \\
      ${args} \\
      --sample ${sample_id} \\
      --seq ${input_reads_arg} \\
      --threads ${task.cpus} \\
      --output ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     mykrobe:
      version: \$(echo \$(mykrobe --version 2>&1) | sed 's/^.*mykrobe v// ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

    stub:
    output = "${sample_id}_mykrobe.csv"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     mykrobe:
      version: \$(echo \$(mykrobe --version 2>&1) | sed 's/^.*mykrobe v// ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
