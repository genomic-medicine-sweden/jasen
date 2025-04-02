process chewbbaca_create_batch_list {
  scratch params.scratch

  input:
    path masked_assembly

  output:
    path "batch_input.list", emit: list

  when:
    task.ext.when

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
