process fastqc {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}*_fastqc/summary.txt")      , emit: summary
    tuple val(sample_id), path("${sample_id}*_fastqc/fastqc_data.txt")  , emit: output
    tuple val(sample_id), path("*.zip")                                 , emit: zip
    tuple val(sample_id), path("*.html")                                , emit: html
    tuple val("${task.process}"), val('fastqc'), eval("echo \$(fastqc --version 2>&1) | sed -r 's/^.*FastQC v//'"), topic: versions, emit: versions_fastqc

    when:
    task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def old_new_pairs = reads instanceof Path || reads.size() == 1 ? [[ reads, "${sample_id}.${reads.extension}" ]] : reads.withIndex().collect { entry, index -> [ entry, "${sample_id}_${index + 1}.${entry.extension}" ] }
    def rename_to     = old_new_pairs*.join(' ').join(' ')
    def renamed_files = old_new_pairs.collect{ _old_name, new_name -> new_name }.join(' ')
    def memory_in_mb  = task.memory ? task.memory.toUnit('MB').toFloat() / task.cpus : null
    def fastqc_memory = memory_in_mb > 10000 ? 10000 : (memory_in_mb < 100 ? 100 : memory_in_mb) // FastQC memory value allowed range (100 - 10000)
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
    """

    stub:
    """
    mkdir ${sample_id}_fastqc

    touch ${sample_id}_fastqc/summary.txt
    touch ${sample_id}_fastqc/fastqc_data.txt
    touch ${sample_id}.zip
    touch ${sample_id}.html
    """
}
