process create_analysis_result {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    path runInfo
    //paths
    tuple val(sampleName), val(quast), val(mlst), val(cgmlst), val(resistance), val(resfinderMeta), val(virulence), val(virulencefinderMeta), val(bracken)

  output:
    path(output)

  script:
    output = "${sampleName}_result.json"
    quastArgs = quast ? "--quast ${quast}" : "" 
    brackenArgs = bracken ? "--kraken ${bracken}" : "" 
    mlstArgs = mlst ? "--mlst ${mlst}" : "" 
    cgmlstArgs = cgmlst ? "--cgmlst ${cgmlst}" : "" 
    resfinderArgs = resistance ? "--resistance ${resistance}" : "" 
    resfinderArgs = resfinderMeta ? "${resfinderArgs} --process-metadata ${resfinderMeta}" : resfinderArgs
    virulenceArgs = virulence ? "--virulence ${virulence}" : "" 
    virulenceArgs = virulencefinderMeta ? "${virulenceArgs} --process-metadata ${virulencefinderMeta}" : virulenceArgs
    """
    prp create-output \\
      --sample-id ${sampleName} \\
      --run-metadata ${runInfo} \\
      ${quastArgs} \\
      ${brackenArgs} \\
      ${mlstArgs} \\
      ${cgmlstArgs} \\
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
