process bracken {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(report)
    path database

  output:
    tuple val(sampleName), path(output)      , emit: output
    tuple val(sampleName), path(outputReport), emit: report
    path "*versions.yml"                     , emit: versions

  script:
    def args = task.ext.args ?: ''
    output = "${sampleName}_bracken.out"
    outputReport = "${sampleName}_bracken.report"
    """
    bracken \\
    ${args} \\
    -d ${database} \\
    -i ${report} \\
    -o ${output} \\
    -w ${outputReport}

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     bracken:
      version: 2.8
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleName}_bracken.out"
    outputReport = "${sampleName}_bracken.report"
    """
    touch $output
    touch $outputReport

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     bracken:
      version: 2.8
      container: ${task.container}
    END_VERSIONS
    """
}
