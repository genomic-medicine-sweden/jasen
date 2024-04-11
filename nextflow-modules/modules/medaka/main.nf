process medaka {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(reads), val(platform) 
    // path assembly 
    tuple val(sampleName), path(assembly)
  output:
    tuple val(sampleName), path("${sampleName}_consensus.fasta"), emit: fasta
    path "*versions.yml"               , emit: versions

  when:
    // task.ext.when && platform == "nanopore"
    platform == "nanopore"

  script:
    def args = task.ext.args ?: ''
    outputDir = params.publishDir ? params.publishDir : 'medaka'
    """
    medaka_consensus -i ${reads} -d ${assembly} -o ${outputDir} ${args}
    mv ${outputDir}/consensus.fasta ${sampleName}_consensus.fasta

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     medaka:
      version: \$(echo \$(medaka --version 2>&1) | sed 's/medaka //')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleName}_consensus.fasta"
    """
    touch $output

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     medaka:
      version: \$(echo \$(medaka --version 2>&1) | sed 's/medaka //')
      container: ${task.container}
    END_VERSIONS
    """
}
