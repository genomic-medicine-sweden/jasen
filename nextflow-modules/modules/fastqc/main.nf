process fastqc {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(reads)

  output:
    tuple val(sampleID), path(summary_output) , emit: summary
    tuple val(sampleID), path(output)         , emit: output
    tuple val(sampleID), path(html_output)    , emit: html
    path "*versions.yml"                      , emit: versions

  script:
    def memory_in_mb = task.memory ? task.memory.toUnit('MB').toFloat() / task.cpus : null
    def fastqc_memory = memory_in_mb > 10000 ? 10000 : (memory_in_mb < 100 ? 100 : memory_in_mb) // FastQC memory value allowed range (100 - 10000)
    def reads_arg = reads.size() == 2 ? "${reads[0]} ${reads[1]}" : "${reads[0]}"
    summary_output = "${sampleID}_fastqc_summary.txt"
    output = "${sampleID}_fastqc_data.txt"
    html_output = "${sampleID}_fastqc.html"
    """
    fastqc \\
      ${args} \\
      --extract \\
      --threads ${task.cpus} \\
      --memory ${fastqc_memory} \\
      -o . \\
      ${reads_arg}

    mv */summary.txt ${summary_output}
    mv */fastqc_data.txt ${output}
    mv *_fastqc.html ${html_output}

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     fastqc:
      version: \$(echo \$(fastqc --version 2>&1) | sed -r 's/^.*FastQC v//')
      container: ${task.container}
    END_VERSIONS
    """

 stub:
    summary_output = "${sampleID}_fastqc_summary.txt"
    output = "${sampleID}_fastqc_data.txt"
    html_output = "${sampleID}_fastqc.html"
    """
    touch ${summary_output}
    touch ${output}
    touch ${html_output}

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     fastqc:
      version: \$(echo \$(fastqc --version 2>&1) | sed -r 's/^.*FastQC v//')
      container: ${task.container}
    END_VERSIONS
    """
}
