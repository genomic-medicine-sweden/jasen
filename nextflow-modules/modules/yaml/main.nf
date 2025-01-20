process create_yaml {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(prp), path(signature), path(ska_index)
    val species

  output:
    tuple val(sample_id), path(output), emit: yaml

  script:
    def prp_dir = task.ext.args_prp ?: ""
    def sourmash_dir = task.ext.args_sourmash ?: ""
    def ska_dir = task.ext.args_ska ?: ""
    output = "${sample_id}_bonsai.yaml"
    """
    #!/usr/bin/env python
    import yaml

    # Define the data
    data = {
        "group_id": "${species}",
        "prp_result": "${prp_dir}/${prp}",
        "minhash_signature": "${sourmash_dir}/${signature}",
        "ska_index": "${ska_dir}/${ska_index}"
    }

    # Write the data to the YAML file
    with open("${output}", "w") as fout:
        yaml.dump(data, fout, default_flow_style=False)
    """

  stub:
    output = "${sample_id}_bonsai.yaml"
    """
    touch ${output}
    """
}
