process flye {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path(output), emit: fasta
    path "*versions.yml"              , emit: versions

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

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     flye:
      version: \$(echo \$(flye --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """

    stub:
    output = "${sample_id}_flye.fasta"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     flye:
      version: \$(echo \$(flye --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """
}
