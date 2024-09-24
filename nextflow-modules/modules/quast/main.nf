process quast {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(assembly)
    path reference

  output:
    tuple val(sampleID), path(output), emit: qc
    path "*versions.yml"             , emit: versions

  when:
    workflow.profile != "streptococcus"

  script:
    def args = task.ext.args ?: ''
    output = "${sampleID}_quast.tsv"
    reference_command = reference ? "-r ${reference}" : ''
    outputDir = "quast_outdir"
    """
    quast.py $args $assembly $reference_command -o $outputDir -t ${task.cpus}
    cp $outputDir/transposed_report.tsv $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     quast:
      version: \$(echo \$(quast.py --version 2>&1) | sed 's/^.*QUAST v//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleID}_quast.tsv"
    """
    touch $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     quast:
      version: \$(echo \$(quast.py --version 2>&1) | sed 's/^.*QUAST v//')
      container: ${task.container}
    END_VERSIONS
    """
}
