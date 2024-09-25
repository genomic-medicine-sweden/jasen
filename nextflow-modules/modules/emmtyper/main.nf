process emmtyper {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(fasta)

  output:
    tuple val(sampleID), path(output), emit: txt
    path "*versions.yml"             , emit: versions

  when:
    workflow.profile == "streptococcus" || workflow.profile == "streptococcus_pyogenes"

  script:
    def args = task.ext.args ?: ''
    output = "${sampleID}_emmtyper.txt"
    """
    emmtyper ${args} --outout ${output} ${fasta}

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     emmtyper:
      version: \$(echo \$(emmtyper --version 2>&1) | sed -r 's/^.*emmtyper // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleID}_emmtyper.txt"
    """
    touch $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     emmtyper:
      version: \$(echo \$(emmtyper --version 2>&1) | sed -r 's/^.*emmtyper // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
