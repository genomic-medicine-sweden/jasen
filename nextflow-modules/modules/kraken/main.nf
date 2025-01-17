process kraken {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(reads)
    path database

  output:
    tuple val(sample_id), path(output), emit: output
    tuple val(sample_id), path(report), emit: report
    path "*versions.yml"              , emit: versions

  when:
    task.ext.when

  script:
    def args = task.ext.args ?: ''
    def reads_arg = reads.size() == 2 ? "${reads[0]} ${reads[1]}" : "${reads[0]}"
    output = "${sample_id}_kraken.out"
    report = "${sample_id}_kraken.report"
    """
    kraken2 \\
    ${args} \\
    --threads ${task.cpus} \\
    --db ${database} \\
    --output ${output} \\
    --report ${report} \\
    ${reads_arg}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     kraken2:
      version: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sample_id}_kraken.out"
    report = "${sample_id}_kraken.report"
    """
    touch ${output}
    touch ${report}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     kraken2:
      version: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
