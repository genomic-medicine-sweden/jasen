process spades {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path(output), emit: fasta
    tuple val("${task.process}"), val('spades'), eval("echo \$(spades.py --version 2>&1) | sed 's/^.*SPAdes genome assembler v//'"), topic: versions, emit: versions_spades

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_reads_arg = reads.size() == 2 ? "-1 ${reads[0]} -2 ${reads[1]}" : "-s ${reads[0]}"
    output_dir = "spades_outdir"
    output = "${sample_id}_spades.fasta"
    """
    spades.py ${args} ${input_reads_arg} -t ${task.cpus} -o ${output_dir}
    mv ${output_dir}/contigs.fasta ${output}
    """

    stub:
    output = "${sample_id}_spades.fasta"
    """
    touch ${output}
    """
}
