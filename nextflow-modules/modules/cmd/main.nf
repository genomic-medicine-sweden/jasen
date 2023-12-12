process export_to_cdm {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(qc), val(sequencingRun), val(limsID)
    val species

  output:
    tuple val(sampleName), path(output), emit: cdm

  script:
    output = "${sampleName}.cdm"
    sequencingRun = sequencingRun ? "-sequencingRun ${sequencingRun}" : ""
    limsID = limsID ? "-lims-id ${limsID}" : ""
    """
    echo ${sequencingRun} \\
         -sample-id ${sampleName} \\
         -assay ${species} \\
         -qc ${qc} \\
         ${limsID} > ${output}
    """

  stub:
    output = "${sampleName}.cdm"
    """
    touch $output
    """
}
