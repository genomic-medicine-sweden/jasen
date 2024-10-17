process kraken {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(reads)
    path database

  output:
    tuple val(sampleID), path(output), emit: output
    tuple val(sampleID), path(report), emit: report
    path "*versions.yml"             , emit: versions

  when:
    params.useKraken

  script:
    def args = task.ext.args ?: ''
    output = "${sampleID}_kraken.out"
    report = "${sampleID}_kraken.report"
    """
    kraken2 \\
    ${args} \\
    --threads ${task.cpus} \\
    --db ${database} \\
    --output ${output} \\
    --report ${report} \\
    ${reads.join(' ')}

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     kraken2:
      version: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleID}_kraken.out"
    report = "${sampleID}_kraken.report"
    """
    touch ${output}
    touch ${report}

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     kraken2:
      version: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
