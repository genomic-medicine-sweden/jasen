process bwa_index {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(fasta), val(platform)

  output:
    tuple val(sample_id), path("${fasta}.*"), emit: index
    path "*versions.yml"                    , emit: versions

  when:
    task.ext.when && platform == "illumina"

  script:
    """
    bwa index ${fasta} ${fasta.baseName}/${fasta}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     bwa:
      version: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    """
    touch ${fasta}.amb
    touch ${fasta}.ann
    touch ${fasta}.bwt
    touch ${fasta}.pac
    touch ${fasta}.sa

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
    tuple val(sample_id), path(reads), val(platform)
    path referenceIdx

  output:
    tuple val(sample_id), path(output), emit: bam
    path "*versions.yml"              , emit: versions

  when:
    task.ext.when && platform == "illumina"

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
