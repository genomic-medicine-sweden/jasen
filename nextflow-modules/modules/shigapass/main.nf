process shigapass {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(assembly)
    path shigapass_db

  output:
    tuple val(sample_id), path(output), emit: csv
    path "*versions.yml"              , emit: versions

  when:
    task.ext.when

  script:
    shigapass_db_arg = shigapass_db ? "-p ${shigapass_db}" : "-p /usr/local/share/shigapass-1.5.0/db/"
    output = "${sample_id}_shigapass.csv"
    """
    echo ${assembly} > batch_input.txt

    ShigaPass.sh \\
    -l batch_input.txt \\
    -o . \\
    ${shigapass_db_arg} \\
    -t ${task.cpus}

    cp ShigaPass_summary.csv ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     shigapass:
      version: \$(echo \$(ShigaPass.sh -v 2>&1) | sed 's/^.*ShigaPass version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

 stub:
    output = "${sample_id}_shigapass.csv"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     shigapass:
      version: \$(echo \$(ShigaPass.sh -v 2>&1) | sed 's/^.*ShigaPass version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
