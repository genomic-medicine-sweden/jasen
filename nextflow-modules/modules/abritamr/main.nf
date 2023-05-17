process abritamr {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(contigs)
    val species
    path database

  output:
    tuple val(sampleName), path(matches)  , optional: true, emit: matches
    tuple val(sampleName), path(partials) , optional: true, emit: partials
    tuple val(sampleName), path(virulence), optional: true, emit: virulence
    path "*versions.yml"                  , emit: versions

  script:
    def mode = task.ext.mode ?: 'run'
    def database_command = database ? "--amrfinder_db ${database}" : ""
    def matches_command = task.ext.mode == 'report' ? "--matches ${matches}" : ""
    def partials_command = task.ext.mode == 'report' ? "--partials ${partials}" : ""
    def args = task.ext.args ?: ''
    matches = "${sampleName}_matches.tsv"
    partials = "${sampleName}_partials.tsv"
    virulence = "${sampleName}_virulence.tsv"
    output = "${sampleName}_amr.out"
    """
    abritamr $mode \\
    --contigs $contigs \\
    --species $species \\
    $database_command \\
    $matches_command \\
    $partials_command \\
    $args

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     abritamr:
      version: \$(echo \$(abritamr --version 2>&1) | sed 's/^.*abritamr // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    """
    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     abritamr:
      version: \$(echo \$(abritamr --version 2>&1) | sed 's/^.*abritamr // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
