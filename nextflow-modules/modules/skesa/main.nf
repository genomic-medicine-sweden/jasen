process skesa {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(reads), val(platform)

  output:
    tuple val(sample_id), path(output), emit: fasta
    path "*versions.yml"              , emit: versions

  when:
    task.ext.when && platform == "illumina"

  script:
    def args = task.ext.args ?: ''
    def input_reads_arg = reads.size() == 2 ? "${reads[0]},${reads[1]}" : "${reads[0]}"
    output = "${sample_id}_skesa.fasta"
    """
    skesa --reads ${input_reads_arg} ${args} > ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     skesa:
      version: \$(echo \$(skesa --version 2>&1) | sed 's/^.*SKESA // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sample_id}_skesa.fasta"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     skesa:
      version: \$(echo \$(skesa --version 2>&1) | sed 's/^.*SKESA // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
