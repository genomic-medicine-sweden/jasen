process ska_build {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path(output), emit: skf
    path "*versions.yml"              , emit: versions

    script:
    def args = task.ext.args ?: ''
    def input_reads_arg = reads.size() == 2 ? "${sample_id}\t${reads[0]}\t${reads[1]}" : "${sample_id}\t${reads[0]}"
    output_basename = "${sample_id}_ska_index"
    output = "${output_basename}.skf"
    """
    echo ${input_reads_arg} > ${sample_id}_input.txt
    ska build ${args} --threads ${task.cpus} -o ${output_basename} -f ${sample_id}_input.txt

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     ska2:
      version: \$(echo \$(ska --version 2>&1) | sed 's/^.*ska // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

    stub:
    output = "${sample_id}_ska_index.skf"
    """
    mkdir ${sample_id}
    touch $output

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     ska2:
      version: \$(echo \$(ska --version 2>&1) | sed 's/^.*ska // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
