process bwa_index {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(reference)

  output:
    tuple val(sampleName), path("${reference}.*"), emit: idx
    path "*versions.yml"                         , emit: versions

  script:
    """
    bwa index ${reference} ${reference.baseName}/${reference}

    cat <<-END_VERSIONS > ${task.process}_versions.yml
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

    cat <<-END_VERSIONS > ${task.process}_versions.yml
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
    tuple val(sampleName), path("${sampleName}.sam") , emit: sam
    path "*versions.yml"                             , emit: versions

  script:
    def args = task.ext.args ?: ''
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    bwa mem \\
      ${args} \\
      -t ${task.cpus} \\
      \${INDEX} \\
      ${reads.join(' ')} > ${sampleName}.sam 

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     bwa:
      version: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    """
    touch ${sampleName}.sam

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     bwa:
      version: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
      container: ${task.container}
    END_VERSIONS
    """
}
