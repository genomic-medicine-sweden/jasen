process quast {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(assembly)
    path reference

  output:
    tuple val(sampleName), path("${assembly.simpleName}.quast.tsv"), emit: qc
    path "*versions.yml"                                           , emit: versions

  script:
    def args = task.ext.args ?: ''
    sampleName = "${assembly.simpleName}"
    reference = reference ? "-r ${reference}" : ''
    outputDir = 'quast_outdir'
    """
    quast.py ${args} ${assembly} ${reference} -o ${outputDir} -t ${task.cpus}
    cp ${outputDir}/transposed_report.tsv ${sampleName}.quast.tsv

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     quast:
      version: \$(echo \$(quast.py --version 2>&1) | sed 's/^.*QUAST v//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    """
    touch ${assembly.simpleName}.quast.tsv

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     quast:
      version: \$(echo \$(quast.py --version 2>&1) | sed 's/^.*QUAST v//')
      container: ${task.container}
    END_VERSIONS
    """
}
