process shigatyper {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.tsv")     , emit: tsv
    tuple val(sample_id), path("${sample_id}-hits.tsv"), optional: true, emit: hits
    path "*versions.yml"                               , emit: versions

    when:
    task.ext.when

    script:
    def is_paired = (reads instanceof List) && reads.size() == 2
    def reads_arg
    if (params.platform == "nanopore") {
        def single_input = (reads instanceof List) ? reads[0] : reads
        reads_arg = "--SE ${single_input} --ont"
    } else if (is_paired) {
        reads_arg = "--R1 ${reads[0]} --R2 ${reads[1]}"
    } else {
        def single_input = (reads instanceof List) ? reads[0] : reads
        reads_arg = "--SE ${single_input}"
    }
    """
    shigatyper \\
        ${reads_arg} \\
        --name ${sample_id}

    cat <<END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     shigatyper:
      version: \$(echo \$(shigatyper --version 2>&1) | sed -n 's/.*ShigaTyper v\\. \\([^,]*\\),.*/\\1/p')
      container: ${task.container}
END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.tsv

    cat <<END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     shigatyper:
      version: \$(echo \$(shigatyper --version 2>&1) | sed -n 's/.*ShigaTyper v\\. \\([^,]*\\),.*/\\1/p')
      container: ${task.container}
END_VERSIONS
    """
}
