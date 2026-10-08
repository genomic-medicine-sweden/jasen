process medaka {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads), path(assembly)

    output:
    tuple val(sample_id), path(output), emit: fasta
    tuple val("${task.process}"), val('medaka'), eval("medaka --version 2>&1 | sed 's/medaka //'"), topic: versions, emit: versions_medaka

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    output_dir = "medaka_outdir"
    output = "${sample_id}_medaka.fasta"
    """
    medaka_consensus -i ${reads} -d ${assembly} -o medaka_tmp ${args}

    medaka_consensus -i ${reads} -d medaka_tmp/consensus.fasta -o ${output_dir} ${args}
    mv ${output_dir}/consensus.fasta ${output}
    """

    stub:
    output = "${sample_id}_medaka.fasta"
    """
    touch ${output}
    """
}
