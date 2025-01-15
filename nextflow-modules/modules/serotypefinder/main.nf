process serotypefinder {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(assembly)
    val databases
    path serotype_db

  output:
    tuple val(sampleID), path(output)     , emit: json
    tuple val(sampleID), path(meta_output), emit: meta
    path "*versions.yml"                  , emit: versions

  when:
    !(workflow.profile in ["mycobacterium_tuberculosis", "staphylococcus_aureus", "streptococcus", "streptococcus_pyogenes"])

  script:
    databases_arg = databases ? "--databases ${databases.join(',')}" : ""
    output = "${sampleID}_serotypefinder.json"
    meta_output = "${sampleID}_serotypefinder_meta.json"
    """
    # Get db version
    DB_HASH=\$(git -C ${serotype_db} rev-parse HEAD)
    JSON_FMT='{"name": "%s", "version": "%s", "type": "%s"}'
    printf "\$JSON_FMT" "serotypefinder" "\$DB_HASH" "database" > ${meta_output}

    # Run serotypefinder
    serotypefinder           \\
    --infile ${assembly}     \\
    ${databases_arg}         \\
    --databasePath ${serotype_db}
    cp data.json ${output}

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     serotypefinder_db:
      version: \$(echo \$DB_HASH | tr -d '\n')
      container: ${task.container}
    END_VERSIONS
    """

 stub:
    output = "${sampleID}_serotypefinder.json"
    meta_output = "${sampleID}_serotypefinder_meta.json"
    """
    DB_HASH=\$(git -C ${serotype_db} rev-parse HEAD)
    touch ${output}
    touch ${meta_output}

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     serotypefinder_db:
      version: \$(echo \$DB_HASH | tr -d '\n')
      container: ${task.container}
    END_VERSIONS
    """
}
