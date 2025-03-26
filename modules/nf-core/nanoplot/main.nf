process nanoplot {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(reads), val(platform) 

  output:
    tuple val(sample_id), path(output), emit: html
    path "*versions.yml"              , emit: versions

  when:
    platform == "nanopore"

  script:
    def args = task.ext.args ?: ''
    output = "${sample_id}_NanoPlot-report.html"
    """
    NanoPlot ${args} --threads ${task.cpus} --prefix ${sample_id}_ --fastq ${reads}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     nanoplot:
      version: \$(echo \$(NanoPlot --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sample_id}_NanoPlot-report.html"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     nanoplot:
      version: \$(echo \$(NanoPlot --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """
}
