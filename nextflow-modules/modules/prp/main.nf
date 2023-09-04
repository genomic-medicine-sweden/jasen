process create_analysis_result {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), val(quast), val(postalignqc), val(mlst), val(cgmlst), val(amr), val(resistance), val(resfinderMeta), val(virulence), val(virulencefinderMeta), val(runInfo), val(mykrobe), val(snippy), val(tbprofiler), val(bracken)

  output:
    path(output)

  script:
    output = "${sampleName}_result.json"
    amrfinderArgs = amr ? "--amr ${amr}" : ""
    brackenArgs = bracken ? "--kraken ${bracken}" : ""
    cgmlstArgs = cgmlst ? "--cgmlst ${cgmlst}" : ""
    mlstArgs = mlst ? "--mlst ${mlst}" : ""
    mykrobeArgs = mykrobe ? "--mykrobe ${mykrobe}" : ""
    postalignqcArgs = postalignqc ? "--quality ${postalignqc}" : "" 
    quastArgs = quast ? "--quast ${quast}" : ""
    resfinderArgs = resistance ? "--resistance ${resistance}" : ""
    resfinderArgs = resfinderMeta ? "${resfinderArgs} --process-metadata ${resfinderMeta}" : resfinderArgs
    runInfoArgs = runInfo ? "--run-metadata ${runInfo}" : ""
    snippyArgs = snippy ? "--snippy ${snippy}" : ""
    tbprofilerArgs = tbprofiler ? "--tbprofiler ${tbprofiler}" : ""
    virulenceArgs = virulence ? "--virulence ${virulence}" : ""
    virulenceArgs = virulencefinderMeta ? "${virulenceArgs} --process-metadata ${virulencefinderMeta}" : virulenceArgs
    """
    prp create-output \\
      --sample-id ${sampleName} \\
      ${amrfinderArgs} \\
      ${brackenArgs} \\
      ${cgmlstArgs} \\
      ${mlstArgs} \\
      ${mykrobeArgs} \\
      ${postalignqcArgs} \\
      ${quastArgs} \\
      ${resfinderArgs} \\
      ${runInfoArgs} \\
      ${snippyArgs} \\
      ${tbprofilerArgs} \\
      ${virulenceArgs} \\
      ${output}
    """

  stub:
    output = "${sampleName}_result.json"
    """
    touch $output
    """
}
