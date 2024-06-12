process skesa {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(reads), val(platform) 

  output:
    tuple val(sampleID), path(output), emit: fasta
    path "*versions.yml"             , emit: versions

  when:
    task.ext.when && platform == "illumina"

  script:
    def args = task.ext.args ?: ''
    def inputData = reads.size() == 2 ? "${reads[0]},${reads[1]}" : "${reads[0]}"
    output = "${sampleID}_skesa.fasta"
    """
    skesa --reads ${inputData} ${args} > ${output}

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     skesa:
      version: \$(echo \$(skesa --version 2>&1) | sed 's/^.*SKESA // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleID}_skesa.fasta"
    """
    touch $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     skesa:
      version: \$(echo \$(skesa --version 2>&1) | sed 's/^.*SKESA // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
