process export_to_cgviz {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    path runInfo
    file meta
    //paths
    tuple val(sampleName), val(quast), val(mlst), val(cgmlst), val(virulence), val(resistance)

  output:
    path(output)

  script:
    output = "${sampleName}_cgviz.json"
    //--kraken ${bracken} \\
    quastArgs = quast ? "--quast ${quast}" : "" 
    mlstArgs = mlst ? "--mlst ${mlst}" : "" 
    cgmlstArgs = cgmlst ? "--cgmlst ${cgmlst}" : "" 
    resfinderArgs = resistance ? "--resistance ${resistance}" : "" 
    virulenceArgs = virulence ? "--virulence ${virulence}" : "" 
    metaArgs = meta ? "--process-metadata  ${meta[1..-1].join(' --process-metadata ')}" : ""
    """
    combine_results.py \\
      --sample-id ${sampleName} \\
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
    output = "${sampleName}_cgviz.json"
    """
    touch $output
    """
}