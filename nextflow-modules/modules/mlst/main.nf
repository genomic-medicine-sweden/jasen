def getAbbrevSpeciesName(fullName) {
  "Convert the full name to the abbreviated version"
  names = fullName.split(' ')
  return names[0][0] + names[1]
}

process mlst {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(assembly)
    val species
    path blastDb

  output:
    tuple val(sampleName), path('*.tsv')  , optional: true, emit: tsv
    tuple val(sampleName), path('*.json') , optional: true, emit: json
    tuple val(sampleName), path('*.novel'), optional: true, emit: novel
    path "*versions.yml"                  , emit: versions

  script:
    def args = task.ext.args ?: ''
    outputName = "${sampleName}.mlst"
    abbrevName = getAbbrevSpeciesName(species)
    blastDbPath = blastDb ? "--blastdb ${blastDb}/mlst.fa" : ""
    //pubmlstDataDir = pubmlstData ? "--datadir ${pubmlstData}" : ""
    //${pubmlstDataDir} \\
    """
    mlst \\
      ${args} \\
      ${blastDbPath} \\
      --scheme  ${abbrevName} \\
      --json ${outputName}.json \\
      --novel ${outputName}.novel \\
      --threads ${task.cpus} \\
      ${assembly}

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     mlst:
      version: \$(echo \$(mlst --version 2>&1) | sed 's/^.*mlst //')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    """
    touch ${sampleName}.tsv
    touch ${sampleName}.json
    touch ${sampleName}.novel

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     mlst:
      version: \$(echo \$(mlst --version 2>&1) | sed 's/^.*mlst //')
      container: ${task.container}
    END_VERSIONS
    """
}
