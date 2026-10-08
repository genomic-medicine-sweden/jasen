process freebayes {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(assembly), path(bam), path(bai) 

    output:
    tuple val(sample_id), path(output), emit: vcf
    tuple val("${task.process}"), val('freebayes'), eval("freebayes --version 2>&1 | sed -r 's/^.*version:[[:space:]]+v// ; s/ .*//'"), topic: versions, emit: versions_freebayes

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    output = "${sample_id}_freebayes.vcf"
    """
    freebayes ${args} -f ${assembly} ${bam} > ${output}
    """

    stub:
    output = "${sample_id}_freebayes.vcf"
    """
    touch ${output}
    """
}
