process minimap2_align {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads)
    path referenceGenomeMmi

    output:
    tuple val(sample_id), path(output), emit: sam
    tuple val("${task.process}"), val('minimap2'), eval("echo \$(minimap2 --version 2>&1)"), topic: versions, emit: versions_minimap2

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    def process = task.process.tokenize(':')[-1]
    output = "${sample_id}_${process}.sam"
    """
    minimap2 ${args} ${referenceGenomeMmi} ${reads} > ${output}
    """

    stub:
    def process = task.process.tokenize(':')[-1]
    output = "${sample_id}_${process}.sam"
    """
    touch "${output}"
    """
}

process minimap2_index {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("*.mmi"), emit: index
    tuple val("${task.process}"), val('minimap2'), eval("echo \$(minimap2 --version 2>&1)"), topic: versions, emit: versions_minimap2

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    minimap2 \\
        -t ${task.cpus} \\
        -d ${fasta.baseName}.mmi \\
        ${args} \\
        ${fasta}
    """

    stub:
    """
    touch ${fasta.baseName}.mmi
    """
}
