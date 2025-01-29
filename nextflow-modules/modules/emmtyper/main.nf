process emmtyper {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(assembly)

  output:
    tuple val(sample_id), path(output), emit: tsv
    path "*versions.yml"              , emit: versions

  when:
    params.species in ["streptococcus", "streptococcus pyogenes"]

  script:
    def args = task.ext.args ?: ''
    output = "${sample_id}_emmtyper.tsv"
    """
    emmtyper ${args} --output ${output} ${assembly}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     emmtyper:
      version: \$(echo \$(emmtyper --version 2>&1) | sed -r 's/^.*emmtyper // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sample_id}_emmtyper.tsv"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     emmtyper:
      version: \$(echo \$(emmtyper --version 2>&1) | sed -r 's/^.*emmtyper // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
