process mask_polymorph_assembly {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(assembly), path(polymorph)

  output:
    tuple val(sampleID), path(output), emit: fasta

  script:
    output = "${sampleID}_mask.fasta"
    """
    error_corr_assembly.pl ${assembly} ${polymorph} > ${output}
    """

  stub:
    output = "${sampleID}_mask.fasta"
    """
    touch $output
    """
}