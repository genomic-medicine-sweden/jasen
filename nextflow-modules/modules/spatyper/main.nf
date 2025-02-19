process spatyper {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(assembly)

  output:
    tuple val(sample_id), path(output), emit: tsv 
    path "*versions.yml"              , emit: versions

  when:
    task.ext.when

  script:
    def args = task.ext.args ?: ''
    output = "${sample_id}_spatyper.tsv"
    """
    spaTyper -f ${assembly} --output ${output} ${args}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     spatyper:
      version: \$(echo \$(spaTyper --version 2>&1) | sed 's/spaTyper //p')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sample_id}_spatyper.tsv"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     spatyper:
      version: \$(echo \$(spatyper --version 2>&1) | sed 's/spaTyper //p')
      container: ${task.container}
    END_VERSIONS
    """
}
