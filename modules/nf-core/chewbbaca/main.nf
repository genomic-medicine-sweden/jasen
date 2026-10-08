process chewbbaca_allelecall {
    tag "${workflow.runName}"
    scratch params.scratch

    input:
    path batch_input
    path schema_dir
    path training_file

    output:
    path('output_dir/results_alleles.tsv'), emit: calls
    tuple val("${task.process}"), val('chewbbaca'), eval("chewie --version 2>/dev/null | sed 's/^.*chewBBACA version: //'"), topic: versions, emit: versions_chewbbaca

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    training_file_arg = training_file ? "--ptf ${training_file}" : "" 
    """
    chewie AlleleCall \\
    -i ${batch_input} \\
    ${args} \\
    --cpu ${task.cpus} \\
    --output-directory output_dir \\
    ${training_file_arg} \\
    --schema-directory ${schema_dir}
    """

    stub:
    """
    mkdir output_dir
    touch output_dir/results_alleles.tsv
    """
}
