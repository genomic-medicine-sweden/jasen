process mykrobe {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(reads)

  output:
    tuple val(sampleName), path(output), emit: json
    path "*versions.yml"               , emit: versions

  script:
    def args = task.ext.args ?: ''
    def inputData = reads.size() == 2 ? "${reads.join(' ')}" : "${reads[0]}"
    output = "${sampleName}_mykrobe.json"
    """
    mykrobe predict \\
      ${args} \\
      --sample ${sampleName} \\
      --seq ${inputData} \\
      --threads ${task.cpus} \\
      --output ${output}

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     mykrobe:
      version: \$(echo \$(mykrobe --version 2>&1) | sed 's/^.*mykrobe v// ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleName}_mykrobe.json"
    """
    touch $output

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     mykrobe:
      version: \$(echo \$(mykrobe --version 2>&1) | sed 's/^.*mykrobe v// ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
