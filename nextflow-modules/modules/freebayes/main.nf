process freebayes {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(fasta)
    tuple path(bam), path(bai) 

  output:
    tuple val(sampleName), path(output), emit: vcf
    path "*versions.yml"               , emit: versions

  script:
    def args = task.ext.args ?: ''
    output = "${sampleName}.vcf"
    """
    freebayes ${args} -f ${fasta} ${bam} > ${output}

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     freebayes:
      version: \$(echo \$(freebayes --version 2>&1) | sed 's/^.*version:  v// ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleName}.vcf"
    """
    touch $output

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     freebayes:
      version: \$(echo \$(freebayes --version 2>&1) | sed 's/^.*version:  v// ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
