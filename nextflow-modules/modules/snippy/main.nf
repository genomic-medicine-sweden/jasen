process snippy {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(reads)
    path reference

  output:
    tuple val(sampleName), path("${sampleName}/snps.vcf"), emit: vcf
    tuple val(sampleName), path("${sampleName}/snps.csv"), emit: csv
    path "*versions.yml"                                 , emit: versions

  when:
    task.ext.when && workflow.profile == "mycobacterium_tuberculosis"

  script:
    def args = task.ext.args ?: ''
    def inputData = reads.size() == 2 ? "--R1 ${reads[0]} --R2 ${reads[1]}" : "--R1 ${reads[0]}"
    """
    snippy \\
      ${args} \\
      ${inputData} \\
      --ref ${reference} \\
      --cpus ${task.cpus} \\
      --ram ${task.memory} \\
      --output ${sampleName}

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     snippy:
      version: \$(echo \$(snippy --version 2>&1) | sed 's/^.*snippy // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    """
    mkdir ${sampleName}
    touch "${sampleName}/snps.{vcf,bed,gff,csv,tab,html,bam,txt}"

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     snippy:
      version: \$(echo \$(snippy --version 2>&1) | sed 's/^.*snippy // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
