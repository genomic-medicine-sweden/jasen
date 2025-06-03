process serotypefinder {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(assembly)
    val databases
    path serotypefinder_db

    output:
    tuple val(sample_id), path(output)     , emit: json
    tuple val(sample_id), path(meta_output), emit: meta
    path "*versions.yml"                   , emit: versions

    when:
    task.ext.when

    script:
    databases_arg = databases ? "--databases ${databases.join(',')}" : ""
    output = "${sample_id}_serotypefinder.json"
    meta_output = "${sample_id}_serotypefinder_meta.json"
    """
    # Get db version
    DB_VERSION=\$(tr -d '\r\n' < ${serotypefinder_db}/VERSION)
    JSON_FMT='{"name": "%s", "version": "%s", "type": "%s"}'
    printf "\$JSON_FMT" "serotypefinder" "\$DB_VERSION" "database" > ${meta_output}

    # Run serotypefinder
    serotypefinder           \\
    --infile ${assembly}     \\
    ${databases_arg}         \\
    --databasePath ${serotypefinder_db}
    cp data.json ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     serotypefinder_db:
      version: \$(echo \$DB_VERSION)
      container: ${task.container}
    END_VERSIONS
    """

 stub:
    output = "${sample_id}_serotypefinder.json"
    meta_output = "${sample_id}_serotypefinder_meta.json"
    """
    DB_VERSION=\$(tr -d '\r\n' < ${serotypefinder_db}/VERSION)
    touch ${output}
    touch ${meta_output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     serotypefinder_db:
      version: \$(echo \$DB_VERSION)
      container: ${task.container}
    END_VERSIONS
    """
}
