process flye {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path(output), emit: fasta
    tuple val("${task.process}"), val('flye'), eval("flye --version 2>&1"), topic: versions, emit: versions_flye

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    def seqmethod = task.ext.seqmethod ?: ''
    def input_reads_arg = reads.size() == 2 ? "${reads[0]} ${reads[1]}" : "${reads[0]}"
    output_dir = "flye_outdir"
    output = "${sample_id}_flye.fasta"
    """
    flye \\
        ${seqmethod} \\
        ${input_reads_arg} \\
        ${args} \\
        --threads ${task.cpus} \\
        --out-dir ${output_dir}

    mv ${output_dir}/assembly.fasta ${output}
    """

    stub:
    output = "${sample_id}_flye.fasta"
    """
    touch ${output}
    """
}
