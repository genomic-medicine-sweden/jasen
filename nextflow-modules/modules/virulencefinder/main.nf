process virulencefinder {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(reads)
    val databases
    path virulenceDb

  output:
    tuple val(sampleName), path(outputFile), emit: json
    tuple val(sampleName), path(metaFile)  , emit: meta
    
  script:
    databasesArgs = databases ? "--databases ${databases.join(',')}" : ""
    outputFile = "${sampleName}_virulencefinder.json"
    metaFile = "${sampleName}_virulencefinder_meta.json"
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
    """

 stub:
    outputFile = "${sampleName}_virulencefinder.json"
    metaFile = "${sampleName}_virulencefinder_meta.json"
    """
    touch $outputFile
    touch $metaFile
    """
}
