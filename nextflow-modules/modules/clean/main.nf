process assembly_trim_clean {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(reads), val(platform)

  output:
    tuple val(sampleID), path(output)

  when:
    platform == "iontorrent"

  script:
    def args = task.ext.args ?: ''
    output = "${sampleID}_cleaned.fastq.gz"
    """
    run_assembly_trimClean.pl --numcpus ${task.cpus} ${args} -i ${reads} -o ${output}
    """

  stub:
    output = "${sampleID}_cleaned.fastq.gz"
    """
    touch $output
    """
}
