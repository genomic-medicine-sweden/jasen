def getSpeciesTaxonName(fullName) {
  "Convert the full name to the abbreviated version"
  names = fullName.split(' ')
  if (fullName == "escherichia coli") {
      return names[0].capitalize()
  }
  return names[0].capitalize() + "_" + names[1]
}

process amrfinderplus {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(assembly)
    val species
    path database

  output:
    tuple val(sampleName), path(output), emit: output
    path "*versions.yml"               , emit: versions

  script:
    def args = task.ext.args ?: ''
    def database_command = database ? "--database ${database}" : ""
    taxonName = getSpeciesTaxonName(species)
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
