process mykrobe {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path(output), emit: csv
    tuple val("${task.process}"), val('mykrobe'), eval("echo \$(mykrobe --version 2>&1) | sed 's/^.*mykrobe v// ; s/ .*//'"), topic: versions, emit: versions_mykrobe

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_reads_arg = reads.size() == 2 ? "${reads.join(' ')}" : "${reads[0]}"
    output = "${sample_id}_mykrobe.csv"
    """
    mykrobe predict \\
      ${args} \\
      --sample ${sample_id} \\
      --seq ${input_reads_arg} \\
      --threads ${task.cpus} \\
      --output ${output}
    """

    stub:
    output = "${sample_id}_mykrobe.csv"
    """
    touch ${output}
    """
}
