process flye {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(reads), val(platform), val(mode) 

  output:
    tuple val(sampleName), path("${sampleName}_assembly.fasta"), emit: fasta
    path "*versions.yml"                                       , emit: versions

  when:
    task.ext.when && platform == "nanopore"

  script:
    def args = task.ext.args ?: ''
    outputDir = params.publishDir ? params.publishDir : 'flye'
    """
    flye ${mode} ${reads} ${args} -o ${outputDir}
    mv ${outputDir}/assembly.fasta ${sampleName}_assembly.fasta

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     flye:
      version: \$(echo \$(flye --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleName}_assembly.fasta"
    """
    touch $output

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     flye:
      version: \$(echo \$(flye --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """
}
