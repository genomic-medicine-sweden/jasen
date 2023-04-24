process kraken {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(reads)
    path database

  output:
    tuple val(sampleName), path(output), emit: output
    tuple val(sampleName), path(report), emit: report
    path "*versions.yml"                    , emit: versions

  script:
    def args = task.ext.args ?: ''
    output = "${sampleName}_kraken.out"
    report = "${sampleName}_kraken.report"
    """
    kraken2 \\
    ${args} \\
    --threads ${task.cpus} \\
    --db ${database} \\
    --output ${output} \\
    --report ${report} \\
    ${reads.join(' ')}

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     kraken2:
      version: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleName}_kraken.out"
    report = "${sampleName}_kraken.report"
    """
    touch ${output}
    touch ${report}

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     kraken2:
      version: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
