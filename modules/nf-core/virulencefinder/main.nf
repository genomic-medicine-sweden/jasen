process virulencefinder {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads)
    val databases
    path virulencefinder_db

    output:
    tuple val(sample_id), path(output)     , emit: json
    tuple val(sample_id), path(meta_output), emit: meta
    path "*versions.yml"                   , emit: versions

    when:
    task.ext.when

    script:
    databases_arg = databases ? "--databases ${databases.join(',')}" : ""
    output = "${sample_id}_virulencefinder.json"
    meta_output = "${sample_id}_virulencefinder_meta.json"
    """
    # Get db version
    DB_VERSION=\$(tr -d '\r\n' < ${virulencefinder_db}/VERSION)
    JSON_FMT='{"name": "%s", "version": "%s", "type": "%s"}'
    printf "\$JSON_FMT" "virulencefinder" "\$DB_VERSION" "database" > ${meta_output}

    # Run virulencefinder
    virulencefinder.py              \\
    --infile ${reads.join(' ')}     \\
    ${databases_arg}                \\
    --databasePath ${virulencefinder_db}
    cp data.json ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     virulencefinder_db:
      version: \$(echo \$DB_VERSION)
      container: ${task.container}
    END_VERSIONS
    """

 stub:
    output = "${sample_id}_virulencefinder.json"
    meta_output = "${sample_id}_virulencefinder_meta.json"
    """
    DB_VERSION=\$(tr -d '\r\n' < ${virulencefinder_db}/VERSION)
    touch ${output}
    touch ${meta_output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     virulencefinder_db:
      version: \$(echo \$DB_VERSION)
      container: ${task.container}
    END_VERSIONS
    """
}
