process export_to_cgviz {
  tag "${sample_id}"
  scratch params.scratch

  input:
    path run_info
    file meta
    tuple val(sample_id), val(quast), val(mlst), val(cgmlst), val(virulence), val(resistance)

  output:
    path(output), emit: json

  script:
    output = "${sample_id}_cgviz.json"
    //--kraken ${bracken} \\
    quast_arg = quast ? "--quast ${quast}" : "" 
    mlst_arg = mlst ? "--mlst ${mlst}" : "" 
    cgmlst_arg = cgmlst ? "--cgmlst ${cgmlst}" : "" 
    resfinder_arg = resistance ? "--resistance ${resistance}" : "" 
    virulence_arg = virulence ? "--virulence ${virulence}" : "" 
    meta_arg = meta ? "--process-metadata  ${meta[1..-1].join(' --process-metadata ')}" : ""
    """
    combine_results.py \\
      --sample-id ${sample_id} \\
      --run-metadata ${run_info} \\
      ${meta_arg} \\
      ${quast_arg} \\
      ${mlst_arg} \\
      ${cgmlst_arg} \\
      ${virulence_arg} \\
      ${resfinder_arg} \\
      ${output}
    """

  stub:
    output = "${sample_id}_cgviz.json"
    """
    touch ${output}
    """
}