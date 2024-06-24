process sourmash {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(assembly)

  output:
    tuple val(sampleID), path(output), emit: signature
    path "*versions.yml"             , emit: versions

  script:
    def args = task.ext.args ?: ''
    output = "${sampleID}.sig"
    """
    sourmash sketch dna $args $assembly -o $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     sourmash:
      version: \$(echo \$(sourmash --version 2>&1) | sed 's/^.*sourmash // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleID}.sig"
    """
    touch $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     sourmash:
      version: \$(echo \$(sourmash --version 2>&1) | sed 's/^.*sourmash // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
