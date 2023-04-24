process freebayes {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(fasta)
    tuple path(bam), path(bai) 

  output:
    tuple val(sampleName), path(output), emit: vcf

  script:
    def args = task.ext.args ?: ''
    output = "${sampleName}.vcf"
    """
    freebayes ${args} -f ${fasta} ${bam} > ${output}
    """

  stub:
    output = "${sampleName}.vcf"
    """
    touch $output
    """
}
