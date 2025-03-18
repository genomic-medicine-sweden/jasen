process minimap2_to_ref {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(reads), val(platform)
    path referenceGenomeMmi
  output:
    tuple val(sample_id), path(output), emit: bam
    path "*versions.yml"              , emit: versions

  when:
    platform == "nanopore"

  script:
    def args = task.ext.args ?: ''
    output = "${sample_id}_minimap2.bam"
    """
    minimap2 ${args} ${referenceGenomeMmi} ${reads} > ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     minimap2:
      version: \$(echo \$(minimap2 --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    """
    touch "${sample_id}_minimap2.bam"

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     minimap2:
      version: \$(echo \$(minimap2 --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """
}

