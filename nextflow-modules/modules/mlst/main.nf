process mlst {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(assembly)
    val scheme
    path pubMlstDb
    path mlstBlastDb

  output:
    tuple val(sampleID), path('*.tsv')  , optional: true, emit: tsv
    tuple val(sampleID), path('*.json') , optional: true, emit: json
    tuple val(sampleID), path('*.novel'), optional: true, emit: novel
    path "*versions.yml"                , emit: versions

  script:
    def args = task.ext.args ?: ''
    outputName = "${sampleID}_mlst"
    schemeArgs = scheme ? "--scheme ${scheme}" : "" 
    pubMlstDbArgs = pubMlstDb ? "--datadir ${pubMlstDb}" : ""
    mlstBlastDbPath = mlstBlastDb ? "--blastdb ${mlstBlastDb}/mlst.fa" : ""
    """
    mlst \\
      ${args} \\
      ${pubMlstDbArgs} \\
      ${schemeArgs} \\
      --json ${outputName}.json \\
      --novel ${outputName}.novel \\
      --threads ${task.cpus} \\
      ${assembly}

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     mlst:
      version: \$(echo \$(mlst --version 2>&1) | sed 's/^.*mlst //')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    outputName = "${sampleID}_mlst"
    """
    touch ${outputName}.tsv
    touch ${outputName}.json
    touch ${outputName}.novel

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     mlst:
      version: \$(echo \$(mlst --version 2>&1) | sed 's/^.*mlst //')
      container: ${task.container}
    END_VERSIONS
    """
}
