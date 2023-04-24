process samtools_view {
  tag "$input"
  scratch params.scratch

  input:
    path input
    path fasta

  output:
    path('*.bam'), optional: true, emit: bam
    path('*.cram'), optional: true, emit: cram
    path "*versions.yml"

  script:
    def reference = fasta ? "--reference ${fasta} -C" : ""
    def prefix = input.simpleName
    def fileType = input.getExtension()
    """
    samtools view $reference ${input} > ${prefix}.${fileType}

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    """
    touch ${sampleName}.bam
    touch ${sampleName}.cram

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}

process samtools_sort {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(input)
    path fasta

  output:
    tuple val(sampleName), path('*.bam') , optional: true, emit: bam
    tuple val(sampleName), path('*.cram'), optional: true, emit: cram
    path "*versions.yml"                 , emit: versions

  script:
    def reference = fasta ? "--reference ${fasta} -O cram" : "-O bam"
    def prefix = input.simpleName
    def fileType = fasta ? "cram" : "bam"
    """
    samtools sort ${reference} -@ $task.cpus -o ${prefix}.sorted.${fileType} ${input}

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    """
    touch ${sampleName}.bam
    touch ${sampleName}.cram

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}

process samtools_index {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(input)

  output:
    tuple val(sampleName), path(output), emit: bai
    path "*versions.yml"               , emit: versions

  script:
    output = "${input}.bai"
    """
    samtools index -@ $task.cpus ${input}

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${input}.bai"
    """
    touch $output

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}

process samtools_faidx {
  tag "$input"
  scratch params.scratch

  input:
    path input

  output:
    path output         , emit: fai
    path "*versions.yml", emit: versions

  script:
    output = "${input}.fai"
    """
    samtools faidx ${input}

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${input}.fai"
    """
    touch $output

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
