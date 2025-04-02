process fastqc {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads), val(platform)

    output:
    tuple val(sample_id), path(summary_output), emit: summary
    tuple val(sample_id), path(output)        , emit: output
    tuple val(sample_id), path(html_output)   , emit: html
    path "*versions.yml"                      , emit: versions

    when:
    platform == "illumina"

    script:
    def args          = task.ext.args ?: ''
    def old_new_pairs = reads instanceof Path || reads.size() == 1 ? [[ reads, "${sample_id}.${reads.extension}" ]] : reads.withIndex().collect { entry, index -> [ entry, "${sample_id}_${index + 1}.${entry.extension}" ] }
    def rename_to     = old_new_pairs*.join(' ').join(' ')
    def renamed_files = old_new_pairs.collect{ _old_name, new_name -> new_name }.join(' ')
    def memory_in_mb  = task.memory ? task.memory.toUnit('MB').toFloat() / task.cpus : null
    def fastqc_memory = memory_in_mb > 10000 ? 10000 : (memory_in_mb < 100 ? 100 : memory_in_mb) // FastQC memory value allowed range (100 - 10000)
    summary_output    = "${sample_id}_fastqc_summary.txt"
    output            = "${sample_id}_fastqc_data.txt"
    html_output       = "${sample_id}_fastqc.html"
    """
    printf "%s %s\\n" ${rename_to} | while read old_name new_name; do
        [ -f "\${new_name}" ] || ln -s \$old_name \$new_name
    done

    fastqc \\
      ${args} \\
      --extract \\
      --threads ${task.cpus} \\
      --memory ${fastqc_memory} \\
      --outdir . \\
      ${renamed_files}

    cp ${sample_id}_1_fastqc/summary.txt ${summary_output}
    cp ${sample_id}_1_fastqc/fastqc_data.txt ${output}
    cp ${sample_id}_1_fastqc.html ${html_output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     fastqc:
      version: \$(echo \$(fastqc --version 2>&1) | sed -r 's/^.*FastQC v//')
      container: ${task.container}
    END_VERSIONS
    """

 stub:
    summary_output = "${sample_id}_fastqc_summary.txt"
    output = "${sample_id}_fastqc_data.txt"
    html_output = "${sample_id}_fastqc.html"
    """
    touch ${summary_output}
    touch ${output}
    touch ${html_output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     fastqc:
      version: \$(echo \$(fastqc --version 2>&1) | sed -r 's/^.*FastQC v//')
      container: ${task.container}
    END_VERSIONS
    """
}
