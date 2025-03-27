process minimap2_align {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(reads), val(platform)
    path referenceGenomeMmi

  output:
    tuple val(sample_id), path(output), emit: sam
    path "*versions.yml"              , emit: versions

  when:
    platform == "nanopore"

  script:
    def args = task.ext.args ?: ''
    def process = task.process.tokenize(':')[-1]
    output = "${sample_id}_${process}.sam"
    """
    minimap2 ${args} ${referenceGenomeMmi} ${reads} > ${output}
    
    cat <<-END_VERSIONS > ${sample_id}_${process}_versions.yml
    ${task.process}:
     minimap2:
      version: \$(echo \$(minimap2 --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    """
    touch "${output}"

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     minimap2:
      version: \$(echo \$(minimap2 --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """
}

process minimap2_index {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(fasta), val(platform)

  output:
    tuple val(sample_id), path("*.mmi"), emit: index
    path "*versions.yml"               , emit: versions

  when:
    platform == "nanopore"

  script:
    def args = task.ext.args ?: ''
    """
    minimap2 \\
        -t ${task.cpus} \\
        -d ${fasta.baseName}.mmi \\
        ${args} \\
        ${fasta}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     minimap2:
      version: \$(echo \$(minimap2 --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    """
    touch ${fasta.baseName}.mmi

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     minimap2:
      version: \$(echo \$(minimap2 --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """
}
