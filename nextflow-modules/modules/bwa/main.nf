process bwa_index {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(reference)

  output:
    tuple val(sample_id), path("${reference}.*"), emit: idx
    path "*versions.yml"                        , emit: versions

  when:
    task.ext.when

  script:
    """
    bwa index ${reference} ${reference.baseName}/${reference}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     bwa:
      version: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    """
    touch ${reference}.amb
    touch ${reference}.ann
    touch ${reference}.bwt
    touch ${reference}.pac
    touch ${reference}.sa

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     bwa:
      version: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
      container: ${task.container}
    END_VERSIONS
    """
}

process bwa_mem {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(reads)
    path referenceIdx

  output:
    tuple val(sample_id), path(output), emit: bam
    path "*versions.yml"             , emit: versions

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

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     bwa:
      version: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
      container: ${task.container}
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sample_id}_bwa.bam"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     bwa:
      version: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
      container: ${task.container}
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
      container: ${task.container}
    END_VERSIONS
    """
}
