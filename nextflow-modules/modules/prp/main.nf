process create_analysis_result {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), val(quast), val(postalignqc), val(mlst), val(cgmlst), val(amr), val(resistance), val(resfinderMeta), val(virulence), val(virulencefinderMeta), val(runInfo), val(bracken)

  output:
    path(output)

  script:
    output = "${sampleName}_result.json"
    quastArgs = quast ? "--quast ${quast}" : ""
    postalignqcArgs = postalignqc ? "--quality ${postalignqc}" : "" 
    brackenArgs = bracken ? "--kraken ${bracken}" : ""
    mlstArgs = mlst ? "--mlst ${mlst}" : ""
    cgmlstArgs = cgmlst ? "--cgmlst ${cgmlst}" : ""
    amrfinderArgs = amr ? "--amr ${amr}" : ""
    resfinderArgs = resistance ? "--resistance ${resistance}" : ""
    resfinderArgs = resfinderMeta ? "${resfinderArgs} --process-metadata ${resfinderMeta}" : resfinderArgs
    virulenceArgs = virulence ? "--virulence ${virulence}" : ""
    virulenceArgs = virulencefinderMeta ? "${virulenceArgs} --process-metadata ${virulencefinderMeta}" : virulenceArgs
    runInfoArgs = runInfo ? "--run-metadata ${runInfo}" : ""
    """
    prp create-output \\
      --sample-id ${sampleName} \\
      ${runInfoArgs} \\
      ${quastArgs} \\
      ${postalignqcArgs} \\
      ${brackenArgs} \\
      ${mlstArgs} \\
      ${cgmlstArgs} \\
      ${amrfinderArgs} \\
      ${virulenceArgs} \\
      ${resfinderArgs} \\
      ${output}
    """

  stub:
    output = "${sampleName}_result.json"
    """
    touch $output
    """
}
