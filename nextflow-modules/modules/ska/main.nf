process ska_build {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(reads)

  output:
    tuple val(sampleID), path(output), emit: skf
    path "*versions.yml"             , emit: versions

  script:
    def args = task.ext.args ?: ''
    def inputData = reads.size() == 2 ? "${sampleID}\t${reads[0]}\t${reads[1]}" : "${sampleID}\t${reads[0]}"
    output_basename = "${sampleID}_ska_index"
    output = "${output_basename}.skf"
    """
    echo ${inputData} > ${sampleID}_input.txt
    ska build $args -o ${output_basename} -f ${sampleID}_input.txt

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     ska2:
      version: \$(echo \$(ska --version 2>&1) | sed 's/^.*ska // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleID}_ska_index.skf"
    """
    mkdir ${sampleID}
    touch $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     ska2:
      version: \$(echo \$(ska --version 2>&1) | sed 's/^.*ska // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
