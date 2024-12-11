process mykrobe {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(reads)

  output:
    tuple val(sampleID), path(output), emit: csv
    path "*versions.yml"             , emit: versions

  when:
    workflow.profile == "mycobacterium_tuberculosis"

  script:
    def args = task.ext.args ?: ''
    def inputData = reads.size() == 2 ? "${reads.join(' ')}" : "${reads[0]}"
    output = "${sampleID}_mykrobe.csv"
    """
    mykrobe predict \\
      ${args} \\
      --sample ${sampleID} \\
      --seq ${inputData} \\
      --threads ${task.cpus} \\
      --output ${output}

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     mykrobe:
      version: \$(echo \$(mykrobe --version 2>&1) | sed 's/^.*mykrobe v// ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleID}_mykrobe.csv"
    """
    touch $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     mykrobe:
      version: \$(echo \$(mykrobe --version 2>&1) | sed 's/^.*mykrobe v// ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
