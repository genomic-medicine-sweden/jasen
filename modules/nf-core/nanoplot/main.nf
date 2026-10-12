process nanoplot {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path(output_html), emit: html
    tuple val(sample_id), path(output_txt),  emit: txt
    tuple val("${task.process}"), val('nanoplot'), eval("NanoPlot --version 2>/dev/null | sed 's/^.*NanoPlot //'"), topic: versions, emit: versions_nanoplot

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    output_html = "${sample_id}_NanoPlot-report.html"
    output_txt = "${sample_id}_NanoStats.txt"
    """
    NanoPlot ${args} --threads ${task.cpus} --prefix ${sample_id}_ --fastq ${reads}
    """

    stub:
    output_html = "${sample_id}_NanoPlot-report.html"
    output_txt = "${sample_id}_NanoStats.txt"
    """
    touch ${output_html} ${output_txt}
    """
}
