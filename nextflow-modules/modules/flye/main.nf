process flye {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(reads), val(platform) 

  output:
    tuple val(sampleID), path(output), emit: fasta
    path "*versions.yml"             , emit: versions

  when:
    platform == "nanopore"

  script:
    def args = task.ext.args ?: ''
    def seqmethod = task.ext.seqmethod ?: ''
    outputDir = "flye_outdir"
    output = "${sampleID}_flye.fasta"
    """
    flye $seqmethod $reads --out-dir $outputDir $args
    mv $outputDir/assembly.fasta $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     flye:
      version: \$(echo \$(flye --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleID}_flye.fasta"
    """
    touch $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     flye:
      version: \$(echo \$(flye --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """
}
