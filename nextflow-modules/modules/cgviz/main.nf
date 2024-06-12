process export_to_cgviz {
  tag "${sampleID}"
  scratch params.scratch

  input:
    path runInfo
    file meta
    tuple val(sampleID), val(quast), val(mlst), val(cgmlst), val(virulence), val(resistance)

  output:
    path(output)

  script:
    output = "${sampleID}_cgviz.json"
    //--kraken ${bracken} \\
    quastArgs = quast ? "--quast ${quast}" : "" 
    mlstArgs = mlst ? "--mlst ${mlst}" : "" 
    cgmlstArgs = cgmlst ? "--cgmlst ${cgmlst}" : "" 
    resfinderArgs = resistance ? "--resistance ${resistance}" : "" 
    virulenceArgs = virulence ? "--virulence ${virulence}" : "" 
    metaArgs = meta ? "--process-metadata  ${meta[1..-1].join(' --process-metadata ')}" : ""
    """
    combine_results.py \\
      --sample-id ${sampleID} \\
      --run-metadata ${runInfo} \\
      ${metaArgs} \\
      ${quastArgs} \\
      ${mlstArgs} \\
      ${cgmlstArgs} \\
      ${virulenceArgs} \\
      ${resfinderArgs} \\
      ${output}
    """

  stub:
    output = "${sampleID}_cgviz.json"
    """
    touch $output
    """
}