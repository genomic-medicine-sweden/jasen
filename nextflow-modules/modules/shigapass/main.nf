process shigapass {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(fasta)
    path shigapassDb

  output:
    tuple val(sampleID), path(outputFile), emit: csv
    path "*versions.yml"                 , emit: versions

  script:
    shigapassDbArgs = shigapassDb ? "-p ${shigapassDb}" : "-p /usr/local/share/shigapass-1.5.0/db/"
    outputFile = "${sampleID}_shigapass.csv"
    """
    echo ${fasta} > batch_input.txt

    ShigaPass.sh \\
    -l batch_input.txt \\
    -o . \\
    ${shigapassDbArgs} \\
    -t ${task.cpus}

    cp ShigaPass_summary.csv ${outputFile}

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     shigapass:
      version: \$(echo \$(ShigaPass.sh -v 2>&1) | sed 's/^.*ShigaPass version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

 stub:
    outputFile = "${sampleID}_shigapass.csv"
    """
    touch $outputFile

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     shigapass:
      version: \$(echo \$(ShigaPass.sh -v 2>&1) | sed 's/^.*ShigaPass version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
