process serotypefinder {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(reads)
    val databases
    path serotypeDb

  output:
    tuple val(sampleName), path(outputFile), emit: json
    tuple val(sampleName), path(metaFile)  , emit: meta
    path "*versions.yml"                   , emit: versions

  script:
    databasesArgs = databases ? "--databases ${databases.join(',')}" : ""
    outputFile = "${sampleName}_serotypefinder.json"
    metaFile = "${sampleName}_serotypefinder_meta.json"
    """
    # Get db version
    DB_HASH=\$(git -C ${serotypeDb} rev-parse HEAD)
    JSON_FMT='{"name": "%s", "version": "%s", "type": "%s"}'
    printf "\$JSON_FMT" "serotypefinder" "\$DB_HASH" "database" > ${metaFile}

    # Run serotypefinder
    serotypefinder              \\
    --infile ${reads.join(' ')}     \\
    ${databasesArgs}                \\
    --databasePath ${serotypeDb}
    cp data.json ${outputFile}

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     serotypefinder_db:
      version: \$(echo \$DB_HASH | tr -d '\n')
      container: ${task.container}
    END_VERSIONS
    """

 stub:
    outputFile = "${sampleName}_serotypefinder.json"
    metaFile = "${sampleName}_serotypefinder_meta.json"
    """
    DB_HASH=\$(git -C ${serotypeDb} rev-parse HEAD)
    touch $outputFile
    touch $metaFile

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     serotypefinder_db:
      version: \$(echo \$DB_HASH | tr -d '\n')
      container: ${task.container}
    END_VERSIONS
    """
}
