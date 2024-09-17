process ska_build {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(assembly)

  output:
    tuple val(sampleID), path(output), emit: skf
    path "*versions.yml"             , emit: versions

  script:
    def args = task.ext.args ?: ''
    output = "${sampleID}/${sampleID}_ska_index.skf"
    """
    ska build $args -o ${output.baseName} $assembly 

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     ska2:
      version: \$(echo \$(ska --version 2>&1) | sed 's/^.*ska // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleID}/${sampleID}_ska_index.skf"
    """
    mkdir ${sampleID}
    touch $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     ska2:
      version: \$(echo \$(ska --version 2>&1) | sed 's/^.*ska // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}

process ska_align {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(skf)

  output:
    tuple val(sampleID), path(output), emit: aln
    path "*versions.yml"             , emit: versions

  script:
    def args = task.ext.args ?: ''
    output = "${sampleID}/${sampleID}_ska_index.aln"
    """
    ska align $args -o ${output} $skf 

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     ska2:
      version: \$(echo \$(ska --version 2>&1) | sed 's/^.*ska // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleID}/${sampleID}_ska_index.aln"
    """
    mkdir ${sampleID}
    touch $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     ska2:
      version: \$(echo \$(ska --version 2>&1) | sed 's/^.*ska // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
