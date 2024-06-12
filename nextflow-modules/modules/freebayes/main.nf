process freebayes {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(fasta)
    tuple path(bam), path(bai) 

  output:
    tuple val(sampleID), path(output), emit: vcf
    path "*versions.yml"             , emit: versions

  script:
    def args = task.ext.args ?: ''
    output = "${sampleID}_freebayes.vcf"
    """
    freebayes ${args} -f ${fasta} ${bam} > ${output}

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     freebayes:
      version: \$(echo \$(freebayes --version 2>&1) | sed 's/^.*version:  v// ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleID}_freebayes.vcf"
    """
    touch $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     freebayes:
      version: \$(echo \$(freebayes --version 2>&1) | sed 's/^.*version:  v// ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
