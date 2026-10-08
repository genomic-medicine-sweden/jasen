process bwa_index {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("bwa"), emit: index
    tuple val("${task.process}"), val('bwa'), eval("bwa 2>&1 | sed 's/^.*Version: //; s/Contact:.*\$//'"), topic: versions, emit: versions_bwa

    when:
    task.ext.when

    script:
    """
    mkdir bwa
    bwa index -p bwa/${fasta.baseName} ${fasta}
    """

    stub:
    """
    mkdir bwa
    touch bwa/${fasta}.amb
    touch bwa/${fasta}.ann
    touch bwa/${fasta}.bwt
    touch bwa/${fasta}.pac
    touch bwa/${fasta}.sa
    """
}

process bwa_mem {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads)
    path index

    output:
    tuple val(sample_id), path(output), emit: bam
    tuple val("${task.process}"), val('bwa'), eval("bwa 2>&1 | sed 's/^.*Version: //; s/Contact:.*\$//'"), topic: versions, emit: versions_bwa
    tuple val("${task.process}"), val('samtools'), eval("samtools --version 2>&1 | sed 's/^.*samtools //; s/Using.*\$//'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    output = "${sample_id}_bwa.bam"
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    bwa mem \\
        ${args} \\
        -t ${task.cpus} \\
        \$INDEX \\
        ${reads.join(' ')} \\
        | samtools sort ${args2} --threads ${task.cpus} -o ${output} -
    """

    stub:
    output = "${sample_id}_bwa.bam"
    """
    touch ${output}
    """
}
