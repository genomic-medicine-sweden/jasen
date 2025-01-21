process sccmec {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(assembly)

  output:
    tuple val(sampleID), path(output), emit: tsv 
    path "*versions.yml"             , emit: versions

  script:
    def args = task.ext.args ?: ''
    output = "${sampleID}_sccmec.tsv"
    outputDir = "sccmec_outdir"
    """
    sccmec --input $assembly --prefix ${sampleID}_sccmec -o $outputDir $args
    cp $outputDir/${sampleID}_sccmec.tsv $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     sccmec:
      version: \$(echo \$(sccmec --version 2>&1) | tr '\n' ';')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleID}_sccmec.tsv"
    """
    touch $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     sccmec:
      version: \$(echo \$(sccmec --version 2>&1) | tr '\n' ';')
      container: ${task.container}
    END_VERSIONS
    """
}
