process hostile {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(reads)

  output:
    tuple val(sampleID), path("${output_dir}/*.fastq.gz")   , emit: reads
    path "*versions.yml"                                    , emit: versions

  when:
    task.ext.when

  script:
    def args = task.ext.args ?: ''
    def reads_arg = reads.size() == 2 ? "--fastq1 ${reads[0]} --fastq2 ${reads[1]}" : "--fastq1 ${reads[0]}"
    output_dir = "hostile_outdir"
    """
    hostile clean ${args} ${reads_arg} --output ${output_dir}

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     hostile:
      version: \$(echo \$(hostile --version 2>&1) )
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output_dir = "hostile_outdir"
    """
    mkdir ${output_dir}
    touch ${output_dir}/${sampleID}_R1.fastq.gz
    touch ${output_dir}/${sampleID}_R2.fastq.gz

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     hostile:
      version: \$(echo \$(hostile --version 2>&1) )
      container: ${task.container}
    END_VERSIONS
    """
}
