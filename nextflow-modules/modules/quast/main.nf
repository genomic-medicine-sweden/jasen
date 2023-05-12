process quast {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(assembly)
    path reference

  output:
    tuple val(sampleName), path(output), emit: qc
    path "*versions.yml"               , emit: versions

  script:
    def args = task.ext.args ?: ''
    output = "${sampleName}_quast.tsv"
    reference_command = reference ? "-r ${reference}" : ''
    outputDir = 'quast_outdir'
    """
    quast.py $args $assembly $reference_command -o $outputDir -t ${task.cpus}
    cp ${outputDir}/transposed_report.tsv $output

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     quast:
      version: \$(echo \$(quast.py --version 2>&1) | sed 's/^.*QUAST v//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleName}_quast.tsv"
    """
    touch $output

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     quast:
      version: \$(echo \$(quast.py --version 2>&1) | sed 's/^.*QUAST v//')
      container: ${task.container}
    END_VERSIONS
    """
}
