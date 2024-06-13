process spades_iontorrent {
  tag "${sampleID}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleID), path(reads), val(platform)

  output:
    tuple val(sampleID), path(output), emit: fasta
    path "*versions.yml"             , emit: versions

  when:
    task.ext.when && platform == "iontorrent"

  script:
    def args = task.ext.args ?: ''
    def inputData = reads.size() == 2 ? "-1 ${reads[0]} -2 ${reads[1]}" : "-s ${reads[0]}"
    outputDir = "spades_outdir"
    output = "${sampleID}_spades.fasta"
    """
    spades.py ${args} ${inputData} -t ${task.cpus} -o $outputDir
    mv $outputDir/contigs.fasta $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     spades:
      version: \$(echo \$(spades.py --version 2>&1) | sed 's/^.*SPAdes genome assembler v//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleID}_spades.fasta"
    """
    touch $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     spades:
      version: \$(echo \$(spades.py --version 2>&1) | sed 's/^.*SPAdes genome assembler v//')
      container: ${task.container}
    END_VERSIONS
    """
}

process spades_illumina {
  tag "${sampleID}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleID), path(reads), val(platform)

  output:
    tuple val(sampleID), path("${sampleID}.fasta"), emit: fasta
    path "*versions.yml"                              , emit: versions

  when:
    task.ext.when && platform == "illumina"

  script:
    def args = task.ext.args ?: ''
    def inputData = reads.size() == 2 ? "-1 ${reads[0]} -2 ${reads[1]}" : "-s ${reads[0]}"
    outputDir = "spades_outdir"
    output = "${sampleID}_spades.fasta"
    """
    spades.py ${args} ${inputData} -t ${task.cpus} -o $outputDir
    mv $outputDir/contigs.fasta $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     spades:
      version: \$(echo \$(spades.py --version 2>&1) | sed 's/^.*SPAdes genome assembler v//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleID}_spades.fasta"
    """
    touch $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     spades:
      version: \$(echo \$(spades.py --version 2>&1) | sed 's/^.*SPAdes genome assembler v//')
      container: ${task.container}
    END_VERSIONS
    """
}
