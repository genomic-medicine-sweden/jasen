process bracken {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(report)
    path database

  output:
    tuple val(sample_id), path(output)       , emit: output
    tuple val(sample_id), path(output_report), emit: report
    path "*versions.yml"                     , emit: versions

  when:
    task.ext.when

  script:
    def args = task.ext.args ?: ''
    output = "${sample_id}_bracken.out"
    output_report = "${sample_id}_bracken.report"
    """
    bracken \\
    ${args} \\
    -d ${database} \\
    -i ${report} \\
    -o ${output} \\
    -w ${output_report}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     bracken:
      version: 2.8
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sample_id}_bracken.out"
    output_report = "${sample_id}_bracken.report"
    """
    touch ${output}
    touch ${output_report}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     bracken:
      version: 2.8
      container: ${task.container}
    END_VERSIONS
    """
}
