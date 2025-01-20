process abritamr {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(contigs)
    val species
    path database

  output:
    tuple val(sample_id), path(matches)  , optional: true, emit: matches
    tuple val(sample_id), path(partials) , optional: true, emit: partials
    tuple val(sample_id), path(virulence), optional: true, emit: virulence
    path "*versions.yml"                , emit: versions

  script:
    def mode = task.ext.mode ?: 'run'
    def database_arg = database ? "--amrfinder_db ${database}" : ""
    def matches_arg = task.ext.mode == 'report' ? "--matches ${matches}" : ""
    def partials_arg = task.ext.mode == 'report' ? "--partials ${partials}" : ""
    def args = task.ext.args ?: ''
    matches = "${sample_id}_matches.tsv"
    partials = "${sample_id}_partials.tsv"
    virulence = "${sample_id}_virulence.tsv"
    output = "${sample_id}_amr.out"
    """
    abritamr ${mode} \\
    --contigs ${contigs} \\
    --species ${species} \\
    ${database_arg} \\
    ${matches_arg} \\
    ${partials_arg} \\
    ${args}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     abritamr:
      version: \$(echo \$(abritamr --version 2>&1) | sed 's/^.*abritamr // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    """
    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     abritamr:
      version: \$(echo \$(abritamr --version 2>&1) | sed 's/^.*abritamr // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
