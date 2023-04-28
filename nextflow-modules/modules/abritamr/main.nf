process abritamr {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(reads)
    path database

  output:
    tuple val(sampleName), path(output), emit: output
    tuple val(sampleName), path(report), emit: report
    path "*versions.yml"                    , emit: versions

  script:
    def args = task.ext.args ?: ''
    output = "${sampleName}_amr.out"
    """
    abritamr run\\
    ${args} \\
    --output ${output}

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     abritamr:
      version: \$(echo \$(abritamr --version 2>&1) | sed 's/^.*abritamr // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleName}_amr.out"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     abritamr:
      version: \$(echo \$(abritamr --version 2>&1) | sed 's/^.*abritamr // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
