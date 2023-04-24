process sambamba_markdup {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(bam)
    path index

  output:
    path outputBam
    path outputIndex

  script:
    outputBam = "${bam.simpleName}_dedup.bam"
    outputIndex = "${outputBam}.bai"
    """
    sambamba markdup --tempdir tmp -t ${task.cpus} ${params.join(' ')} ${bam} ${outputBam}
    """
}
