process shigatyper {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.tsv")     , emit: tsv
    tuple val(sample_id), path("${sample_id}-hits.tsv"), optional: true, emit: hits
    tuple val("${task.process}"), val('shigatyper'), eval("echo \$(shigatyper --version 2>&1) | sed 's/^.*ShigaTyper //'"), topic: versions, emit: versions_shigatyper

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
    """

    stub:
    """
    touch ${sample_id}.tsv
    """
}
