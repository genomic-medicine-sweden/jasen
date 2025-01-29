process snippy {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(reads)
    path reference

  output:
    tuple val(sample_id), path(output)                 , emit: vcf
    tuple val(sample_id), path("${sample_id}/snps.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}/snps.csv"), emit: csv
    path "*versions.yml"                               , emit: versions

  when:
    workflow.profile == "mycobacterium_tuberculosis"

  script:
    def args = task.ext.args ?: ''
    def input_reads_arg = reads.size() == 2 ? "--R1 ${reads[0]} --R2 ${reads[1]}" : "--R1 ${reads[0]}"
    output = "${sample_id}_snippy.vcf"
    """
    snippy \\
      ${args} \\
      ${input_reads_arg} \\
      --ref ${reference} \\
      --cpus ${task.cpus} \\
      --outdir ${sample_id}

    cp ${sample_id}/snps.vcf ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     snippy:
      version: \$(echo \$(snippy --version 2>&1) | sed 's/^.*snippy // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sample_id}_snippy.vcf"
    """
    mkdir ${sample_id}
    touch ${output}
    touch ${sample_id}/snps.{vcf,bed,gff,csv,tab,html,bam,txt}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     snippy:
      version: \$(echo \$(snippy --version 2>&1) | sed 's/^.*snippy // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
