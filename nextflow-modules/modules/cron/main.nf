process copy_to_cron {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), val(prp)

  output:
    tuple val(sampleName), path(output), emit: json

  when:
    params.cronCopy

  script:
    output = "${sampleName}_result.json"
    """
    cp ${prp} ${output}
    """

  stub:
    output = "${sampleName}_result.json"
    """
    touch $output
    """
}
