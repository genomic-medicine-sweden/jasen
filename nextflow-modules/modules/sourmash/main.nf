process sourmash {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(assembly)

  output:
    tuple val(sampleName), path(output), emit: signature
    path "*versions.yml"               , emit: versions

  when:
    task.ext.when

  script:
    def args = task.ext.args ?: ''
    output = "${sampleName}.sig"
    """
    sourmash sketch dna $args $assembly -o $output

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     sourmash:
      version: \$(echo \$(sourmash --version 2>&1) | sed 's/^.*sourmash // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleName}.sig"
    """
    touch $output

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     sourmash:
      version: \$(echo \$(sourmash --version 2>&1) | sed 's/^.*sourmash // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
