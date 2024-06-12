process chewbbaca_allelecall {
  tag "${workflow.runName}"
  scratch params.scratch

  input:
    val sampleID
    path batchInput
    path schemaDir
    path trainingFile

  output:
    val sampleID                          , emit: sampleID
    path('output_dir/results_alleles.tsv'), emit: calls
    path "*versions.yml"                  , emit: versions

  script:
    def args = task.ext.args ?: ''
    trainingFile = trainingFile ? "--ptf ${trainingFile}" : "" 
    """
    chewie AlleleCall \\
    -i ${batchInput} \\
    ${args} \\
    --cpu ${task.cpus} \\
    --output-directory output_dir \\
    ${trainingFile} \\
    --schema-directory ${schemaDir}

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
    path maskedAssembly

  output:
    path "batch_input.list", emit: list

  script:
    output = "batch_input.list"
    """
    realpath $maskedAssembly > $output
    """

  stub:
    output = "batch_input.list"
    """
    touch $output
    """
}

process chewbbaca_split_results {
  tag "${sampleID}"
  scratch params.scratch

  input:
    each sampleID
    path input

  output:
    tuple val(sampleID), path(output), emit: output

  script:
    output = "${sampleID}_chewbbaca.out"
    """
    head -1 ${input} > ${output}
    grep ${sampleID} ${input} >> ${output}
    """

  stub:
    output = "${sampleID}_chewbbaca.out"
    """
    touch $output
    """
}
