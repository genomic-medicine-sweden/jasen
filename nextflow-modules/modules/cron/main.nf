process copy_to_cron {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(yaml), path(cdm)

  output:
    tuple val(sampleName), path(yamlOutput), emit: yaml
    tuple val(sampleName), path(cdmOutput) , emit: cdm

  when:
    params.cronCopy

  script:
    yamlOutput = "${sampleName}_bonsai.yaml"
    cdmOutput = "${sampleName}.cdm"
    """
    cp ${yaml} ${yamlOutput}
    cp ${cdm} ${cdmOutput}
    """

  stub:
    yamlOutput = "${sampleName}_bonsai.yaml"
    cdmOutput = "${sampleName}.cdm"
    """
    touch $yamlOutput $cdmOutput
    """
}
