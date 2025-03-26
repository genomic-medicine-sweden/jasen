process chewbbaca_split_results {
  tag "${sample_id}"
  scratch params.scratch

  input:
    each sample_id
    path input

  output:
    tuple val(sample_id), path(output), emit: output

  when:
    task.ext.when

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
