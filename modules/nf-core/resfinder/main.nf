include { check_taxon } from '../../../methods/check_taxon.nf'

process resfinder {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads)
    val species
    val resfinder_db
    val pointfinder_db

    output:
    tuple val(sample_id), path(output)     , emit: json
    tuple val(sample_id), path(meta_output), emit: meta
    path output_gene                       , emit: gene_table
    path output_point                      , optional: true, emit: point_table
    path "*versions.yml"                   , emit: versions

    when:
    task.ext.when

    script:
    def resfinder_arg = resfinder_db ? "--acquired --db_path_res ${resfinder_db}" : ""
    def pointfinder_arg = pointfinder_db ? "--point --db_path_point ${pointfinder_db}" : ""
    def species_name = check_taxon(species)
    def species_arg = species_name ? "--species '${species_name}'" : ""
    def nanopore_arg = task.ext.nanopore_args ?: ''
    output = "${sample_id}_resfinder.json"
    meta_output = "${sample_id}_resfinder_meta.json"
    output_gene = "${sample_id}_pheno_table.txt"
    output_point = "${sample_id}_point_table.txt"
    """
    # Get db version
    RES_DB_VERSION=\$(tr -d '\r\n' < ${resfinder_db}/VERSION)
    POINT_DB_VERSION=\$(tr -d '\r\n' < ${pointfinder_db}/VERSION)
    JSON_FMT='[{"name": "%s", "version": "%s", "type": "%s"},{"name": "%s", "version": "%s", "type": "%s"}]'
    printf "\$JSON_FMT" "resfinder" \$RES_DB_VERSION "database" "pointfinder" \$POINT_DB_VERSION "database" > ${meta_output}

    # Run resfinder
    python -m resfinder                           \\
    --inputfastq ${reads.join(' ')}               \\
    ${species_arg}                                \\
    ${resfinder_arg}                              \\
    ${pointfinder_arg}                            \\
    ${nanopore_arg}                               \\
    --out_json std_format_under_development.json  \\
    --outputPath .

    cp std_format_under_development.json ${output}
    cp pheno_table.txt ${output_gene}
    if [ -f 'PointFinder_results.txt' ]; then
      cp PointFinder_results.txt ${output_point}
    fi

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     resfinder:
      version: \$(echo \$(python -m resfinder --version 2>&1) )
      container: ${task.container}
     resfinder_db:
      version: \$(echo \$DB_VERSION)
      container: ${task.container}
     pointfinder_db:
      version: \$(echo \$DB_VERSION)
      container: ${task.container}
    END_VERSIONS
    """

    stub:
    output = "${sample_id}_resfinder.json"
    meta_output = "${sample_id}_resfinder_meta.json"
    output_gene = "${sample_id}_pheno_table.txt"
    output_point = "${sample_id}_point_table.txt"
    """
    RES_DB_VERSION=\$(tr -d '\r\n' < ${resfinder_db}/VERSION)
    POINT_DB_VERSION=\$(tr -d '\r\n' < ${pointfinder_db}/VERSION)
    JSON_FMT='[{"name": "%s", "version": "%s", "type": "%s"},{"name": "%s", "version": "%s", "type": "%s"}]'
    printf "\$JSON_FMT" "resfinder" \$RES_DB_VERSION "database" "pointfinder" \$POINT_DB_VERSION "database" > ${meta_output}
    touch ${output}
    touch ${output_gene}
    touch ${output_point}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     resfinder:
      version: \$(echo \$(python -m resfinder --version 2>&1) )
      container: ${task.container}
     resfinder_db:
      version: \$(echo \$DB_VERSION)
      container: ${task.container}
     pointfinder_db:
      version: \$(echo \$DB_VERSION)
      container: ${task.container}
    END_VERSIONS
    """
}
