process bracken {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(report)
    path database

  output:
    tuple val(sampleName), path(output)      , emit: output
    tuple val(sampleName), path(outputReport), emit: report
    path "*versions.yml"                          , emit: versions

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

    cat <<-END_VERSIONS > ${task.process}_versions.yml
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

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     bracken:
      version: 2.8
      container: ${task.container}
    END_VERSIONS
    """
}
