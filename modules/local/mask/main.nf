process mask_polymorph_assembly {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(assembly), path(polymorph)

  output:
    tuple val(sample_id), path(output), emit: fasta

  when:
    task.ext.when

  script:
    output = "${sample_id}_mask.fasta"
    """
    error_corr_assembly.pl ${assembly} ${polymorph} > ${output}
    """

  stub:
    output = "${sample_id}_mask.fasta"
    """
    touch ${output}
    """
}