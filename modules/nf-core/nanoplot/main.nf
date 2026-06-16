process nanoplot {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path(output_html), emit: html
    tuple val(sample_id), path(output_txt),  emit: txt
    path "*versions.yml",                    emit: versions

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    output_html = "${sample_id}_NanoPlot-report.html"
    output_txt = "${sample_id}_NanoStats.txt"
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
    output_html = "${sample_id}_NanoPlot-report.html"
    output_txt = "${sample_id}_NanoStats.txt"
    """
    touch ${output_html} ${output_txt}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     nanoplot:
      version: \$(echo \$(NanoPlot --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """
}
