process copy_to_cron {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(yaml), path(cdm)

  output:
    tuple val(sampleID), path(yamlOutput), emit: yaml
    tuple val(sampleID), path(cdmOutput) , emit: cdm

  when:
    params.cronCopy

  script:
    yamlOutput = "${sampleID}_bonsai.yaml"
    cdmOutput = "${sampleID}.cdm"
    """
    cp ${yaml} ${yamlOutput}
    cp ${cdm} ${cdmOutput}
    """

  stub:
    yamlOutput = "${sampleID}_bonsai.yaml"
    cdmOutput = "${sampleID}.cdm"
    """
    touch $yamlOutput $cdmOutput
    """
}
