process bwa_index {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(reference)

  output:
    tuple val(sampleName), path("${reference}.*"), emit: idx
    path "*versions.yml"                         , emit: versions

  when:
    workflow.profile != "mycobacterium_tuberculosis"

  script:
    """
    bwa index ${reference} ${reference.baseName}/${reference}

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
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

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     bwa:
      version: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
      container: ${task.container}
    END_VERSIONS
    """
}

process bwa_mem {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(reads)
    path referenceIdx

  output:
    tuple val(sampleName), path(output) , emit: bam
    path "*versions.yml"                , emit: versions

  when:
    workflow.profile != "mycobacterium_tuberculosis"

  script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    output = "${sampleName}_bwa.bam"
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    bwa mem \\
        $args \\
        -t $task.cpus \\
        \$INDEX \\
        ${reads.join(' ')} \\
        | samtools sort $args2 --threads ${task.cpus} -o $output -

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
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
    """
    touch ${sampleName}.bam

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
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
