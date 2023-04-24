process post_align_qc {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(bam)
    path bai
    path reference

  output:
    tuple val(sampleName), path(output), emit: qc

  script:
    output = "${sampleName}_bwa.qc"
    """
    postaln_qc.pl ${bam} ${reference} ${sampleName} ${task.cpus} > ${output}
    """

  stub:
    output = "${sampleName}_bwa.qc"
    """
    touch $output
    """
}
