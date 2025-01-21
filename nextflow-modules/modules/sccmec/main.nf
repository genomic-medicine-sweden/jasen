process sccmec {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(assembly)

  output:
    tuple val(sample_id), path(output), emit: tsv 
    path "*versions.yml"              , emit: versions

  script:
    def args = task.ext.args ?: ''
    output = "${sample_id}_sccmec.tsv"
    outputDir = "sccmec_outdir"
    """
    sccmec --input ${assembly} --prefix ${sample_id}_sccmec -o ${outputDir} ${args}
    cp ${outputDir}/${sample_id}_sccmec.tsv ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     sccmec:
      version: \$(echo \$(sccmec --version 2>&1) | sed -n 's/.*sccmec_targets, version //p')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sample_id}_sccmec.tsv"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     sccmec:
      version: \$(echo \$(sccmec --version 2>&1) | sed -n 's/.*sccmec_targets, version //p')
      container: ${task.container}
    END_VERSIONS
    """
}
