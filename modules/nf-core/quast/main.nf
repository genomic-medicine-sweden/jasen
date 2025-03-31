process quast {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(assembly)
    path reference

  output:
    tuple val(sample_id), path(output), emit: tsv
    path "*versions.yml"              , emit: versions

  script:
    def args = task.ext.args ?: ''
    output = "${sample_id}_quast.tsv"
    reference_command = reference ? "-r ${reference}" : ""
    output_dir = "quast_outdir"
    """
    quast.py ${args} ${assembly} ${reference_command} -o ${output_dir} -t ${task.cpus}
    cp ${output_dir}/transposed_report.tsv ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     quast:
      version: \$(echo \$(quast.py --version 2>&1) | sed 's/^.*QUAST v//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sample_id}_quast.tsv"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     quast:
      version: \$(echo \$(quast.py --version 2>&1) | sed 's/^.*QUAST v//')
      container: ${task.container}
    END_VERSIONS
    """
}
