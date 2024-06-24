process bracken {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(report)
    path database

  output:
    tuple val(sampleID), path(output)      , emit: output
    tuple val(sampleID), path(outputReport), emit: report
    path "*versions.yml"                   , emit: versions

  script:
    def args = task.ext.args ?: ''
    output = "${sampleID}_bracken.out"
    outputReport = "${sampleID}_bracken.report"
    """
    bracken \\
    ${args} \\
    -d ${database} \\
    -i ${report} \\
    -o ${output} \\
    -w ${outputReport}

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     bracken:
      version: 2.8
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleID}_bracken.out"
    outputReport = "${sampleID}_bracken.report"
    """
    touch $output
    touch $outputReport

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     bracken:
      version: 2.8
      container: ${task.container}
    END_VERSIONS
    """
}
