process mask_polymorph_assembly {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(assembly), path(polymorph)

  output:
    tuple val(sampleName), path(output), emit: fasta

  script:
    output = "${sampleName}_masked.fasta"
    """
    error_corr_assembly.pl ${assembly} ${polymorph} > ${output}
    """

  stub:
    output = "${sampleName}_masked.fasta"
    """
    touch $output
    """
}