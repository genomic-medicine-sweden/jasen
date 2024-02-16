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
    JSON_FMT='{"name": "%s", "version": "%s", "type": "%s"}'
    printf "\$JSON_FMT" "serotypefinder" "ada62c62a7fa74032448bb2273d1f7045c59fdda" "database" > ${metaFile}

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
    touch $outputFile
    touch $metaFile

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     serotypefinder_db:
      version: \$(echo ada62c62a7fa74032448bb2273d1f7045c59fdda | tr -d '\n')
      container: ${task.container}
    END_VERSIONS
    """
}
