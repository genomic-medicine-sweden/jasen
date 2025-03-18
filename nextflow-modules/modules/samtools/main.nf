process samtools_view {
  tag "${input}"
  scratch params.scratch

  input:
    path input
    path fasta

  output:
    path('*.bam'), optional: true , emit: bam
    path('*.cram'), optional: true, emit: cram
    path "*versions.yml"          , emit: versions
  
  when:
    task.ext.when

  script:
    def reference_arg = fasta ? "--reference ${fasta} -C" : ""
    def prefix = input.simpleName
    def file_ext = input.getExtension()
    """
    samtools view ${reference_arg} ${input} > ${prefix}.${file_ext}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    """
    touch ${sample_id}.bam
    touch ${sample_id}.cram

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}

process samtools_sort {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(input)

  output:
    tuple val(sample_id), path(output), emit: bam
    path "*versions.yml"              , emit: versions

  script:
    output = "${sample_id}.bam"
    """
    samtools sort -@ ${task.cpus} -o ${output} ${input}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    """
    touch ${sample_id}.bam

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}

process samtools_index {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(input)

  output:
    tuple val(sample_id), path(output), emit: bai
    path "*versions.yml"              , emit: versions

  when:
    task.ext.when

  script:
    output = "${input}.bai"
    """
    samtools index -@ ${task.cpus} ${input}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${input}.bai"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}

process samtools_faidx {
  tag "${input}"
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

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${input}.fai"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}

process samtools_coverage {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(input)

  output:
    tuple val(sample_id), path(output), emit: txt
    path "*versions.yml"              , emit: versions

  script:
    def args = task.ext.args ?: ''
    output = "${sample_id}_mapcoverage.txt"
    output_dir = "postmapqc"
    """
    samtools coverage -o ${sample_id}_mapcoverage.txt ${args} ${input}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
        ${task.process}:
        samtools:
          version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
          container: ${task.container}
        END_VERSIONS
    """

  stub:
    """
    touch "${sample_id}_mapcoverage.txt"

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     samtools:
      version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
      container: ${task.container}
    END_VERSIONS
    """
}