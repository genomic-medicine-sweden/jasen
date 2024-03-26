process create_yaml {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(prp), path(signature)
    val species

  output:
    tuple val(sampleName), path(output), emit: yaml

  script:
    def prpDir = task.ext.args_prp ?: ""
    def sourmashDir = task.ext.args_sourmash ?: ""
    output = "${sampleName}_bonsai.yaml"
    """
    #!/usr/bin/env python
    import yaml

    # Define the data
    data = {
        "group_id": "${species}",
        "prp_result": "${prpDir}/${prp}",
        "minhash_signature": "${sourmashDir}/${signature}"
    }

    # Write the data to the YAML file
    with open("$output", "w") as fout:
        yaml.dump(data, fout, default_flow_style=False)
    """

  stub:
    output = "${sampleName}_bonsai.yaml"
    """
    touch $output
    """
}
