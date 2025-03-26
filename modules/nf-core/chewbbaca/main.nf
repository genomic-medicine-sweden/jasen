process chewbbaca_allelecall {
  tag "${workflow.runName}"
  scratch params.scratch

  input:
    path batch_input
    path schema_dir
    path training_file

  output:
    path('output_dir/results_alleles.tsv'), emit: calls
    path "*versions.yml"                  , emit: versions

  when:
    task.ext.when

  script:
    def args = task.ext.args ?: ''
    training_file = training_file ? "--ptf ${training_file}" : "" 
    """
    chewie AlleleCall \\
    -i ${batch_input} \\
    ${args} \\
    --cpu ${task.cpus} \\
    --output-directory output_dir \\
    ${training_file} \\
    --schema-directory ${schema_dir}

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     chewBBACA:
      version: \$(echo \$(chewie --version 2>&1) | sed 's/^.*chewBBACA version: //')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    """
    mkdir output_dir
    touch output_dir/results_alleles.tsv

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     chewBBACA:
      version: \$(echo \$(chewie --version 2>&1) | sed 's/^.*chewBBACA version: //')
      container: ${task.container}
    END_VERSIONS
    """
}
