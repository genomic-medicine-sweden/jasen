process create_analysis_result {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), val(quast), val(postalignqc), val(mlst), val(cgmlst), val(amr), val(resistance), val(resfinderMeta), val(virulence), val(virulencefinderMeta), val(runInfo), val(mykrobe), val(tbprofiler), val(bracken)

  output:
    tuple val(sampleName), path(output), emit: json

  script:
    output = "${sampleName}_result.json"
    amrfinderArgs = amr ? "--amrfinder ${amr}" : ""
    brackenArgs = bracken ? "--kraken ${bracken}" : ""
    cgmlstArgs = cgmlst ? "--cgmlst ${cgmlst}" : ""
    mlstArgs = mlst ? "--mlst ${mlst}" : ""
    mykrobeArgs = mykrobe ? "--mykrobe ${mykrobe}" : ""
    postalignqcArgs = postalignqc ? "--quality ${postalignqc}" : "" 
    quastArgs = quast ? "--quast ${quast}" : ""
    resfinderArgs = resistance ? "--resfinder ${resistance}" : ""
    resfinderArgs = resfinderMeta ? "${resfinderArgs} --process-metadata ${resfinderMeta}" : resfinderArgs
    runInfoArgs = runInfo ? "--run-metadata ${runInfo}" : ""
    tbprofilerArgs = tbprofiler ? "--tbprofiler ${tbprofiler}" : ""
    virulenceArgs = virulence ? "--virulencefinder ${virulence}" : ""
    virulenceArgs = virulencefinderMeta ? "${virulenceArgs} --process-metadata ${virulencefinderMeta}" : virulenceArgs
    """
    prp create-bonsai-input \\
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
      ${tbprofilerArgs} \\
      ${virulenceArgs} \\
      --output ${output}
    """

  stub:
    output = "${sampleName}_result.json"
    """
    touch $output
    """
}

process create_cdm_input {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), val(quast), val(postalignqc), val(cgmlst)

  output:
    tuple val(sampleName), path(output), emit: json

  script:
    output = "${sampleName}_qc_result.json"
    cgmlstArgs = cgmlst ? "--cgmlst ${cgmlst}" : ""
    postalignqcArgs = postalignqc ? "--quality ${postalignqc}" : "" 
    quastArgs = quast ? "--quast ${quast}" : ""
    """
    prp create-cdm-input \\
      ${cgmlstArgs} \\
      ${postalignqcArgs} \\
      ${quastArgs} \\
      --output ${output}
    """

  stub:
    output = "${sampleName}_qc_result.json"
    """
    touch $output
    """
}