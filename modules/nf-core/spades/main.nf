process spades {
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
    def input_reads_arg = reads.size() == 2 ? "-1 ${reads[0]} -2 ${reads[1]}" : "-s ${reads[0]}"
    output_dir = "spades_outdir"
    output = "${sample_id}_spades.fasta"
    """
    spades.py ${args} ${input_reads_arg} -t ${task.cpus} -o ${output_dir}
    mv ${output_dir}/contigs.fasta ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     spades:
      version: \$(echo \$(spades.py --version 2>&1) | sed 's/^.*SPAdes genome assembler v//')
      container: ${task.container}
    END_VERSIONS
    """

    stub:
    output = "${sample_id}_spades.fasta"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     spades:
      version: \$(echo \$(spades.py --version 2>&1) | sed 's/^.*SPAdes genome assembler v//')
      container: ${task.container}
    END_VERSIONS
    """
}
