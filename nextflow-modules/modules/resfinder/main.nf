process resfinder {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(reads)
    val species
    val resfinderDb
    val pointfinderDb

  output:
    tuple val(sampleName), path(outputFileJson), emit: json
    tuple val(sampleName), path(metaFile)      , emit: meta
    path outputFileGene                        , emit: geneTable
    path outputFilePoint                       , emit: pointTable
    path "*versions.yml"                       , emit: versions

  when:
    task.ext.when && workflow.profile != "mycobacterium_tuberculosis"

  script:
    def resfinderFinderParams = pointfinderDb ? "--acquired --db_path_res ${resfinderDb}" : ""
    def pointFinderParams = pointfinderDb ? "--point --db_path_point ${pointfinderDb}" : ""
    def speciesArgs = species ? "--species '${species}'" : ""
    outputFileJson = "${sampleName}_resfinder.json"
    metaFile = "${sampleName}_resfinder_meta.json"
    outputFileGene = "${sampleName}_pheno_table.txt"
    outputFilePoint = "${sampleName}_point_table.txt"
    """
    # Get db version
    RES_DB_VERSION=\$(cat ${resfinderDb}/VERSION | tr -d '\r' | tr -d '\n')
    POINT_DB_VERSION=\$(cat ${pointfinderDb}/VERSION | tr -d '\n')
    JSON_FMT='[{"name": "%s", "version": "%s", "type": "%s"},{"name": "%s", "version": "%s", "type": "%s"}]'
    printf "\$JSON_FMT" "resfinder" \$RES_DB_VERSION "database" "pointfinder" \$POINT_DB_VERSION "database" > $metaFile

    # Run resfinder
    python -m resfinder             \\
    --inputfastq ${reads.join(' ')} \\
    ${speciesArgs}                  \\
    ${resfinderFinderParams}        \\
    ${pointFinderParams}            \\
    --out_json std_format_under_development.json  \\
    --outputPath .

    cp std_format_under_development.json ${outputFileJson}
    cp pheno_table.txt ${outputFileGene}
    cp PointFinder_results.txt ${outputFilePoint}

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     resfinder:
      version: \$(echo \$(python -m resfinder --version 2>&1) )
      container: ${task.container}
     resfinder_db:
      version: \$(cat ${resfinderDb}/VERSION | tr -d '\n')
      container: ${task.container}
     pointfinder_db:
      version: \$(cat ${pointfinderDb}/VERSION | tr -d '\n')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    outputFileJson = "${sampleName}_resfinder.json"
    metaFile = "${sampleName}_resfinder_meta.json"
    outputFileGene = "${sampleName}_pheno_table.txt"
    outputFilePoint = "${sampleName}_point_table.txt"
    """
    RES_DB_VERSION=\$(cat ${resfinderDb}/VERSION | tr -d '\r' | tr -d '\n')
    POINT_DB_VERSION=\$(cat ${pointfinderDb}/VERSION | tr -d '\n')
    JSON_FMT='[{"name": "%s", "version": "%s", "type": "%s"},{"name": "%s", "version": "%s", "type": "%s"}]'
    printf "\$JSON_FMT" "resfinder" \$RES_DB_VERSION "database" "pointfinder" \$POINT_DB_VERSION "database" > $metaFile
    touch $outputFileJson
    touch $outputFileGene
    touch $outputFilePoint

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     resfinder:
      version: \$(echo \$(python -m resfinder --version 2>&1) )
      container: ${task.container}
     resfinder_db:
      version: \$(cat ${resfinderDb}/VERSION | tr -d '\n')
      container: ${task.container}
     pointfinder_db:
      version: \$(cat ${pointfinderDb}/VERSION | tr -d '\n')
      container: ${task.container}
    END_VERSIONS
    """
}
