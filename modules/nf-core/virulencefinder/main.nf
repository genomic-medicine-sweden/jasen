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
    tuple val("${task.process}"), val('virulencefinder_db'), eval("echo \$DB_VERSION"), topic: versions, emit: versions_virulencefinder_db
    tuple val("${task.process}"), val('virulencefinder'), eval("echo \$(python -m virulencefinder --version 2>&1)"), topic: versions, emit: versions_virulencefinder

    when:
    task.ext.when

    script:
    databases_arg = databases ? "--databases ${databases.join(',')}" : ""
    def nanopore_arg = task.ext.nanopore_args ?: ''
    output = "${sample_id}_virulencefinder.json"
    meta_output = "${sample_id}_virulencefinder_meta.json"
    """
    # Get db version
    DB_VERSION=\$(tr -d '\r\n' < ${virulencefinder_db}/VERSION)
    JSON_FMT='{"name": "%s", "version": "%s", "type": "%s"}'
    printf "\$JSON_FMT" "virulencefinder" "\$DB_VERSION" "database" > ${meta_output}

    # Run virulencefinder
    python -m virulencefinder            \\
    --inputfastq ${reads.join(' ')}      \\
    ${databases_arg}                     \\
    ${nanopore_arg}                      \\
    --databasePath ${virulencefinder_db} \\
    --out_json ${output} \\
    --outputPath .
    """

 stub:
    output = "${sample_id}_virulencefinder.json"
    meta_output = "${sample_id}_virulencefinder_meta.json"
    """
    DB_VERSION=\$(tr -d '\r\n' < ${virulencefinder_db}/VERSION)
    touch ${output}
    touch ${meta_output}
    """
}
