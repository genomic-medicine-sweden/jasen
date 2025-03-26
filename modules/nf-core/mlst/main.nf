process mlst {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(assembly)
    val scheme
    path pbmlst_db
    path blast_db

  when:
    task.ext.when

  output:
    tuple val(sample_id), path('*.tsv')  , optional: true, emit: tsv
    tuple val(sample_id), path('*.json') , optional: true, emit: json
    tuple val(sample_id), path('*.novel'), optional: true, emit: novel
    path "*versions.yml"                 , emit: versions

  script:
    def args = task.ext.args ?: ''
    outputName = "${sample_id}_mlst"
    scheme_arg = scheme ? "--scheme ${scheme}" : "" 
    pubmlst_db_arg = pbmlst_db ? "--datadir ${pbmlst_db}" : ""
    blast_db_arg = blast_db ? "--blastdb ${blast_db}/mlst.fa" : ""
    """
    mlst \\
      ${args} \\
      ${pubmlst_db_arg} \\
      ${scheme_arg} \\
      --json ${outputName}.json \\
      --novel ${outputName}.novel \\
      --threads ${task.cpus} \\
      ${assembly}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     mlst:
      version: \$(echo \$(mlst --version 2>&1) | sed 's/^.*mlst //')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    outputName = "${sample_id}_mlst"
    """
    touch ${outputName}.tsv
    touch ${outputName}.json
    touch ${outputName}.novel

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     mlst:
      version: \$(echo \$(mlst --version 2>&1) | sed 's/^.*mlst //')
      container: ${task.container}
    END_VERSIONS
    """
}
