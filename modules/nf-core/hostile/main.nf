process hostile {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads)
    path hostile_dir
    val hostile_idx

    output:
    tuple val(sample_id), path("${output_dir}/*.fastq.gz"), emit: reads
    path "*versions.yml"                                  , emit: versions

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_reads_arg = reads.size() == 2 ? "--fastq1 ${reads[0]} --fastq2 ${reads[1]}" : "--fastq1 ${reads[0]}"
    output_dir = "hostile_outdir"
    """
    export HOSTILE_CACHE_DIR=${hostile_dir}
    mkdir ${output_dir}

    hostile \\
        clean \\
        ${args} \\
        --threads ${task.cpus} \\
        ${input_reads_arg} \\
        --index ${hostile_dir}/${hostile_idx} \\
        --reorder \\
        --airplane \\
        --output ${output_dir}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     hostile:
      version: \$(echo \$(hostile --version 2>&1) )
      container: ${task.container}
    END_VERSIONS
    """

    stub:
    output_dir = "hostile_outdir"
    """
    mkdir ${output_dir}
    touch ${output_dir}/${sample_id}_R1.fastq.gz
    touch ${output_dir}/${sample_id}_R2.fastq.gz

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     hostile:
      version: \$(echo \$(hostile --version 2>&1) )
      container: ${task.container}
    END_VERSIONS
    """
}
