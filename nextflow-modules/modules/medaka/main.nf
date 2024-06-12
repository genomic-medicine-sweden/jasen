process medaka {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(reads), val(platform)
    tuple val(sampleID), path(assembly)

  output:
    tuple val(sampleID), path(output), emit: fasta
    path "*versions.yml"             , emit: versions

  when:
    platform == "nanopore"

  script:
    def args = task.ext.args ?: ''
    outputDir = 'medaka_outdir'
    output = "${sampleID}_medaka.fasta"
    """
    medaka_consensus -i $reads -d $assembly -o medaka_tmp $args

    medaka_consensus -i $reads -d medaka_tmp/consensus.fasta -o $outputDir $args
    mv $outputDir/consensus.fasta $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     medaka:
      version: \$(echo \$(medaka --version 2>&1) | sed 's/medaka //')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleID}_medaka.fasta"
    """
    touch $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     medaka:
      version: \$(echo \$(medaka --version 2>&1) | sed 's/medaka //')
      container: ${task.container}
    END_VERSIONS
    """
}
