process resfinder {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

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
  
  script:
    def resfinderFinderParams = pointfinderDb ? "--acquired --db_path_res ${resfinderDb}" : ""
    def pointFinderParams = pointfinderDb ? "--point --db_path_point ${pointfinderDb}" : ""
    def speciesArgs = species ? "--species '${species}'" : ""
    outputFileJson = "resfinder_${sampleName}.json"
    metaFile = "resfinder_meta_${sampleName}.json"
    outputFileGene = "pheno_table_${sampleName}.txt"
    outputFilePoint = "point_table_${sampleName}.txt"
    """
    # Get db version
    RES_DB_HASH=\$(git -C ${resfinderDb} rev-parse HEAD)
    POINT_DB_HASH=\$(git -C ${pointfinderDb} rev-parse HEAD)
    JSON_FMT='[{"name": "%s", "version": "%s", "type": "%s"},{"name": "%s", "version": "%s", "type": "%s"}]'
    printf "\$JSON_FMT" "resfinder" \$RES_DB_HASH "database" "pointfinder" \$POINT_DB_HASH "database" > $metaFile

    # Run resfinder
    python -m resfinder             \\
    --inputfastq ${reads.join(' ')} \\
    ${speciesArgs}                   \\
    ${resfinderFinderParams}        \\
    ${pointFinderParams}            \\
    --out_json std_format_under_development.json  \\
    --outputPath .

    cp std_format_under_development.json ${outputFileJson}
    cp pheno_table.txt ${outputFileGene}
    cp PointFinder_results.txt ${outputFilePoint}

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     resfinder:
      version: \$(echo \$(python -m resfinder --version 2>&1) )
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    outputFileJson = "resfinder_${sampleName}.json"
    metaFile = "resfinder_meta_${sampleName}.json"
    outputFileGene = "pheno_table_${sampleName}.txt"
    outputFilePoint = "point_table_${sampleName}.txt"
    """
    touch $outputFileJson
    touch $metaFile
    touch $outputFileGene
    touch $outputFilePoint

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     resfinder:
      version: \$(echo \$(python -m resfinder --version 2>&1) )
      container: ${task.container}
    END_VERSIONS
    """
}
