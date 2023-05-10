process amrfinderplus {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(assembly)
    path database

  output:
    tuple val(sampleName), path(output), emit: output
    path "*versions.yml"               , emit: versions

  script:
    def args = task.ext.args ?: ''
    def database_command = database ? "--database ${database}" : ""
    output = "${sampleName}_amr.out"
    """
    amrfinder \\
    --nucleotide $assembly \\
    $database_command \\
    $args \\
    --output $output

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     amrfinderplus:
      version: \$(echo \$(amrfinder --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleName}_amr.out"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     amrfinderplus:
      version: \$(echo \$(amrfinder --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """
}
