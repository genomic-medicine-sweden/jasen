process freebayes {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(assembly)
    tuple val(sample_id), path(bam), path(bai) 

    output:
    tuple val(sample_id), path(output), emit: vcf
    path "*versions.yml"              , emit: versions

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    output = "${sample_id}_freebayes.vcf"
    """
    freebayes ${args} -f ${assembly} ${bam} > ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     freebayes:
      version: \$(echo \$(freebayes --version 2>&1) | sed -r 's/^.*version:\s+v// ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

    stub:
    output = "${sample_id}_freebayes.vcf"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     freebayes:
      version: \$(echo \$(freebayes --version 2>&1) | sed -r 's/^.*version:\s+v// ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
