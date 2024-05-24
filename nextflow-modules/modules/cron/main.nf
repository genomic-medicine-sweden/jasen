process copy_to_cron {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), val(prp), val(cdm)

  output:
    tuple val(sampleName), path(prpOutput), emit: json
    tuple val(sampleName), path(cdmOutput), emit: cdm

  when:
    params.cronCopy

  script:
    prpOutput = "${sampleName}_result.json"
    cdmOutput = "${sampleName}.cdm"
    """
    cp ${prp} ${prpOutput}
    cp ${cdm} ${cdmOutput}
    """

  stub:
    prpOutput = "${sampleName}_result.json"
    cdmOutput = "${sampleName}.cdm"
    """
    touch $prpOutput $cdmOutput
    """
}
