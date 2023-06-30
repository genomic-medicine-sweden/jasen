process amrfinderplus {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(assembly)
    path database

  output:
    tuple val(sampleName), path(output), emit: output
    path "*versions.yml"               , emit: versions

  when:
    task.ext.when && workflow.profile != "mycobacterium_tuberculosis"

  script:
    def args = task.ext.args ?: ''
    def database_command = database ? "--database ${database}" : ""
    output = "${sampleName}_amrfinder.out"
    """
    amrfinder \\
    --nucleotide $assembly \\
    $database_command \\
    $args \\
    --output $output

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     amrfinderplus:
      version: \$(echo \$(amrfinder --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleName}_amrfinder.out"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     amrfinderplus:
      version: \$(echo \$(amrfinder --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """
}
