process medaka {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads), val(platform)
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path(output), emit: fasta
    path "*versions.yml"              , emit: versions

    when:
    platform == "nanopore"

    script:
    def args = task.ext.args ?: ''
    output_dir = "medaka_outdir"
    output = "${sample_id}_medaka.fasta"
    """
    medaka_consensus -i ${reads} -d ${assembly} -o medaka_tmp ${args}

    medaka_consensus -i ${reads} -d medaka_tmp/consensus.fasta -o ${output_dir} ${args}
    mv ${output_dir}/consensus.fasta ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     medaka:
      version: \$(echo \$(medaka --version 2>&1) | sed 's/medaka //')
      container: ${task.container}
    END_VERSIONS
    """

    stub:
    output = "${sample_id}_medaka.fasta"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     medaka:
      version: \$(echo \$(medaka --version 2>&1) | sed 's/medaka //')
      container: ${task.container}
    END_VERSIONS
    """
}
