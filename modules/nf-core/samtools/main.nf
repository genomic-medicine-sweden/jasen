process samtools_view {
    tag "${input}"
    scratch params.scratch

    input:
    path input
    path fasta

    output:
    path('*.bam'), optional: true , emit: bam
    path('*.cram'), optional: true, emit: cram
    tuple val("${task.process}"), val('samtools'), eval("samtools --version 2>&1 | sed 's/^.*samtools // ; s/ .*//'"), topic: versions, emit: versions_samtools
  
    when:
    task.ext.when

    script:
    def reference_arg = fasta ? "--reference ${fasta} -C" : ""
    def prefix = input.simpleName
    def file_ext = input.getExtension()
    """
    samtools view ${reference_arg} ${input} > ${prefix}.${file_ext}
    """

    stub:
    """
    touch ${sample_id}.bam
    touch ${sample_id}.cram
    """
}

process samtools_sort {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(input)

    output:
    tuple val(sample_id), path(output), emit: bam
    tuple val("${task.process}"), val('samtools'), eval("samtools --version 2>&1 | sed 's/^.*samtools // ; s/ .*//'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when

    script:
    output = "${input.baseName}.bam"
    """
    samtools sort -@ ${task.cpus} -O bam -o ${output} ${input}
    """

    stub:
    output = "${input.baseName}.bam"
    """
    touch "${output}"
    """
}

process samtools_index {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(input)

    output:
    tuple val(sample_id), path(output), emit: bai
    tuple val("${task.process}"), val('samtools'), eval("samtools --version 2>&1 | sed 's/^.*samtools // ; s/ .*//'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when

    script:
    output = "${input}.bai"
    """
    samtools index -@ ${task.cpus} ${input}
    """

    stub:
    output = "${input}.bai"
    """
    touch ${output}
    """
}

process samtools_faidx {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path(output), emit: fai
    tuple val("${task.process}"), val('samtools'), eval("samtools --version 2>&1 | sed 's/^.*samtools // ; s/ .*//'"), topic: versions, emit: versions_samtools

    script:
    output = "${fasta}.fai"
    """
    samtools faidx ${fasta}
    """

    stub:
    output = "${fasta}.fai"
    """
    touch ${output}
    """
}

process samtools_coverage {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(input)

    output:
    tuple val(sample_id), path(output), emit: txt
    tuple val("${task.process}"), val('samtools'), eval("samtools --version 2>&1 | sed 's/^.*samtools //; s/Using.*\$//'"), topic: versions, emit: versions_samtools

    script:
    def args = task.ext.args ?: ''
    output = "${input.baseName}_mapcoverage.txt"
    """
    samtools coverage -o ${output} ${args} ${input}
    """

    stub:
    output = "${input.baseName}_mapcoverage.txt"
    """
    touch "${output}"
    """
}

process samtools_stats {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path(output), emit: stats
    tuple val("${task.process}"), val('samtools'), eval("samtools --version 2>&1 | sed 's/^.*samtools //; s/Using.*\$//'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    output = "${sample_id}.stats"
    """
    samtools stats --threads ${task.cpus} ${args} ${bam} > ${output}
    """

    stub:
    output = "${sample_id}.stats"
    """
    touch ${output}
    """
}

process samtools_bedcov {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(bam), path(bai)
    path bed

    output:
    tuple val(sample_id), path(output), emit: coverage
    tuple val("${task.process}"), val('samtools'), eval("samtools --version 2>&1 | sed 's/^.*samtools //; s/Using.*\$//'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    output = "${sample_id}.bedcov.tsv"
    """
    samtools bedcov ${args} ${bed} ${bam} > ${output}
    """

    stub:
    output = "${sample_id}.bedcov.tsv"
    """
    touch ${output}
    """
}
