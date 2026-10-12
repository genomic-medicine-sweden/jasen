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
    tuple val("${task.process}"), val('serotypefinder'), val('2.0.2'), topic: versions, emit: versions_serotypefinder
    tuple val("${task.process}"), val('serotypefinder_db'), eval("echo \$DB_VERSION"), topic: versions, emit: versions_serotypefinder_db

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
    """

 stub:
    output = "${sample_id}_serotypefinder.json"
    meta_output = "${sample_id}_serotypefinder_meta.json"
    """
    DB_VERSION=\$(tr -d '\r\n' < ${serotypefinder_db}/VERSION)
    touch ${output}
    touch ${meta_output}
    """
}
