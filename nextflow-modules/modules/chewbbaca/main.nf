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

process chewbbaca_create_batch_list {
  scratch params.scratch

  input:
    path masked_assembly

  output:
    path "batch_input.list", emit: list

  script:
    output = "batch_input.list"
    """
    realpath ${masked_assembly} > ${output}
    """

  stub:
    output = "batch_input.list"
    """
    touch ${output}
    """
}

process chewbbaca_split_results {
  tag "${sample_id}"
  scratch params.scratch

  input:
    each sample_id
    path input

  output:
    tuple val(sample_id), path(output), emit: output

  script:
    output = "${sample_id}_chewbbaca.out"
    """
    head -1 ${input} > ${output}
    grep ${sample_id} ${input} >> ${output}
    """

  stub:
    output = "${sample_id}_chewbbaca.out"
    """
    touch ${output}
    """
}
