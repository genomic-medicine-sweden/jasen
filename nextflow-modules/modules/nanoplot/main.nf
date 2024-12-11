process nanoplot {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(reads), val(platform) 

  output:
    tuple val(sampleID), path(output), emit: html
    path "*versions.yml"             , emit: versions

  when:
    platform == "nanopore"

  script:
    def args = task.ext.args ?: ''
    output = "${sampleID}_NanoPlot-report.html"
    """
    NanoPlot $args --prefix $sampleID --fastq $reads

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     nanoplot:
      version: \$(echo \$(NanoPlot --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleID}_NanoPlot-report.html"
    """
    touch $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     nanoplot:
      version: \$(echo \$(NanoPlot --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """
}
