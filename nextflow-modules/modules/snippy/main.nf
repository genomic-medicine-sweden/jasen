process snippy {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(reads)
    path reference

  output:
    tuple val(sampleID), path(output)                , emit: vcf
    tuple val(sampleID), path("${sampleID}/snps.bam"), emit: bam
    tuple val(sampleID), path("${sampleID}/snps.csv"), emit: csv
    path "*versions.yml"                             , emit: versions

  script:
    def args = task.ext.args ?: ''
    def inputData = reads.size() == 2 ? "--R1 ${reads[0]} --R2 ${reads[1]}" : "--R1 ${reads[0]}"
    output = "${sampleID}_snippy.vcf"
    """
    snippy \\
      ${args} \\
      ${inputData} \\
      --ref ${reference} \\
      --cpus ${task.cpus} \\
      --outdir ${sampleID}

    cp ${sampleID}/snps.vcf $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     snippy:
      version: \$(echo \$(snippy --version 2>&1) | sed 's/^.*snippy // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleID}_snippy.vcf"
    """
    mkdir ${sampleID}
    touch $output
    touch "${sampleID}/snps.{vcf,bed,gff,csv,tab,html,bam,txt}"

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     snippy:
      version: \$(echo \$(snippy --version 2>&1) | sed 's/^.*snippy // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
