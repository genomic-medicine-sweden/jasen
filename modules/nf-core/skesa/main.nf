process skesa {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path(output), emit: fasta
    tuple val("${task.process}"), val('skesa'), eval("echo \$(skesa --version 2>&1) | sed 's/^.*SKESA // ; s/ .*//'"), topic: versions, emit: versions_skesa

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_reads_arg = reads.size() == 2 ? "${reads[0]},${reads[1]}" : "${reads[0]}"
    output = "${sample_id}_skesa.fasta"
    """
    skesa --cores ${task.cpus} --memory ${task.memory} --reads ${input_reads_arg} ${args} > ${output}
    """

    stub:
    output = "${sample_id}_skesa.fasta"
    """
    touch ${output}
    """
}
