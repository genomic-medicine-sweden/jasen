process assembly_trim_clean {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(reads), val(platform)

  output:
    tuple val(sample_id), path(output), emit: reads

  when:
    platform == "iontorrent"

  script:
    def args = task.ext.args ?: ''
    output = "${sample_id}_cleaned.fastq.gz"
    """
    run_assembly_trimClean.pl --numcpus ${task.cpus} ${args} -i ${reads} -o ${output}
    """

  stub:
    output = "${sample_id}_cleaned.fastq.gz"
    """
    touch ${output}
    """
}
