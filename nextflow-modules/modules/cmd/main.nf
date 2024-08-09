process export_to_cdm {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(qc), val(sequencingRun), val(limsID), val(sampleName)
    val species

  output:
    tuple val(sampleID), path(output), emit: cdm

  script:
    output = "${sampleID}.cdmpy"
    sequencingRun = sequencingRun ? "--sequencing-run ${sequencingRun}" : ""
    limsID = limsID ? "--lims-id ${limsID}" : ""
    """
    echo ${sequencingRun} \\
         --sample-id ${sampleName} \\
         --assay ${species} \\
         --qc ${params.outdir}/${params.speciesDir}/${params.cdmDir}/${qc} \\
         ${limsID} > ${output}
    """

  stub:
    output = "${sampleID}.cdmpy"
    """
    touch $output
    """
}
