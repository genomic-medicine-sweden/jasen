process spades_iontorrent {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(reads), val(platform)

  output:
    tuple val(sampleName), path("${sampleName}.fasta"), emit: fasta
    path "*versions.yml"                              , emit: versions

  when:
    platform == "iontorrent"

  script:
    def args = task.ext.args ?: ''
    def inputData = reads.size() == 2 ? "-1 ${reads[0]} -2 ${reads[1]}" : "-s ${reads[0]}"
    outputDir = params.publishDir ? params.publishDir : 'spades'
    """
    spades.py ${args} ${inputData} -t ${task.cpus} -o ${outputDir}
    mv ${outputDir}/contigs.fasta ${sampleName}.fasta

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     spades:
      version: \$(echo \$(spades.py --version 2>&1) | sed 's/^.*SPAdes genome assembler v//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleName}.fasta"
    """
    touch $output

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     spades:
      version: \$(echo \$(spades.py --version 2>&1) | sed 's/^.*SPAdes genome assembler v//')
      container: ${task.container}
    END_VERSIONS
    """
}

process spades_illumina {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(reads), val(platform)

  output:
    tuple val(sampleName), path("${sampleName}.fasta"), emit: fasta
    path "*versions.yml"                              , emit: versions

  when:
    task.ext.when && platform == "illumina"

  script:
    def args = task.ext.args ?: ''
    def inputData = reads.size() == 2 ? "-1 ${reads[0]} -2 ${reads[1]}" : "-s ${reads[0]}"
    outputDir = params.publishDir ? params.publishDir : 'spades'
    """
    spades.py ${args} ${inputData} -t ${task.cpus} -o ${outputDir}
    mv ${outputDir}/contigs.fasta ${sampleName}.fasta

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     spades:
      version: \$(echo \$(spades.py --version 2>&1) | sed 's/^.*SPAdes genome assembler v//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleName}.fasta"
    """
    touch $output

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     spades:
      version: \$(echo \$(spades.py --version 2>&1) | sed 's/^.*SPAdes genome assembler v//')
      container: ${task.container}
    END_VERSIONS
    """
}
