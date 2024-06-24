process virulencefinder {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(reads)
    val databases
    path virulenceDb

  output:
    tuple val(sampleID), path(outputFile), emit: json
    tuple val(sampleID), path(metaFile)  , emit: meta
    path "*versions.yml"                   , emit: versions

  script:
    databasesArgs = databases ? "--databases ${databases.join(',')}" : ""
    outputFile = "${sampleID}_virulencefinder.json"
    metaFile = "${sampleID}_virulencefinder_meta.json"
    """
    # Get db version
    DB_HASH=\$(git -C ${virulenceDb} rev-parse HEAD)
    JSON_FMT='{"name": "%s", "version": "%s", "type": "%s"}'
    printf "\$JSON_FMT" "virulencefinder" "\$DB_HASH" "database" > ${metaFile}

    # Run virulencefinder
    virulencefinder.py              \\
    --infile ${reads.join(' ')}     \\
    ${databasesArgs}                \\
    --databasePath ${virulenceDb}
    cp data.json ${outputFile}

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     virulencefinder_db:
      version: \$(echo \$DB_HASH | tr -d '\n')
      container: ${task.container}
    END_VERSIONS
    """

 stub:
    outputFile = "${sampleID}_virulencefinder.json"
    metaFile = "${sampleID}_virulencefinder_meta.json"
    """
    DB_HASH=\$(git -C ${virulenceDb} rev-parse HEAD)
    touch $outputFile
    touch $metaFile

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     virulencefinder_db:
      version: \$(echo \$DB_HASH | tr -d '\n')
      container: ${task.container}
    END_VERSIONS
    """
}
