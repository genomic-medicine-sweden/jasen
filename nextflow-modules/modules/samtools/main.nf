process samtools_view {
  tag "$input"
  scratch params.scratch

  input:
    path input
    path fasta

  output:
    path('*.bam'), optional: true , emit: bam
    path('*.cram'), optional: true, emit: cram
    path "*versions.yml"          , emit: versions
  
  when:
    workflow.profile != "mycobacterium_tuberculosis"

  script:
    def reference = fasta ? "--reference ${fasta} -C" : ""
    def prefix = input.simpleName
    def fileType = input.getExtension()
    """
    samtools view $reference ${input} > ${prefix}.${fileType}

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    """
    touch ${sampleID}.bam
    touch ${sampleID}.cram

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}

process samtools_sort {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(input)

  output:
    tuple val(sampleID), path(output), emit: bam
    path "*versions.yml"             , emit: versions

  script:
    output = "${sampleID}.bam"
    """
    samtools sort -@ $task.cpus -o ${output} ${input}

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    """
    touch ${sampleID}.bam

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}

process samtools_index {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(input)

  output:
    tuple val(sampleID), path(output), emit: bai
    path "*versions.yml"             , emit: versions

  script:
    output = "${input}.bai"
    """
    samtools index -@ $task.cpus ${input}

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
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

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
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

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
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

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
