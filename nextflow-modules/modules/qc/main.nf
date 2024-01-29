process post_align_qc {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(bam)
    path bai
    path reference

  output:
    tuple val(sampleName), path(output), emit: qc

  script:
    output = "${sampleName}_bwa.qc"
    """
    prp create-qc-result --bam ${bam} --reference ${reference} --sample-id ${sampleName} --cpus ${task.cpus} --output ${output}
    """

  stub:
    output = "${sampleName}_bwa.qc"
    """
    touch $output
    """
}
