process chewbbaca_allelecall {
  tag "${workflow.runName}"
  scratch params.scratch

  input:
    val sampleName
    path batchInput
    path schemaDir
    path trainingFile

  output:
    val sampleName                        , emit: sampleName
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
    touch /projects/mikro/nobackup/gms/jasen2/chewbbaca_temp.txt
    
    for f in $maskedAssembly
    do 
      echo "Processing $f" > /projects/mikro/nobackup/gms/jasen2/chewbbaca_temp.txt
    done

    
    """
//realpath $maskedAssembly > $output

  //stub:
    //output = "batch_input.list"
    
    //touch $output
    
}

process chewbbaca_split_results {
  tag "${sampleName}"
  scratch params.scratch

  input:
    each sampleName
    path input

  output:
    tuple val(sampleName), path(output), emit: output

  script:
    output = "${sampleName}_chewbbaca.out"
    """
    head -1 ${input} > ${output}
    grep ${sampleName} ${input} >> ${output}
    """

  stub:
    output = "${sampleName}_chewbbaca.out"
    """
    touch $output
    """
}
