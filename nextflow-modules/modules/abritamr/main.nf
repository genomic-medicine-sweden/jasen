process abritamr {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(contigs)
    val species
    path database

  output:
    tuple val(sampleID), path(matches)  , optional: true, emit: matches
    tuple val(sampleID), path(partials) , optional: true, emit: partials
    tuple val(sampleID), path(virulence), optional: true, emit: virulence
    path "*versions.yml"                , emit: versions

  script:
    def mode = task.ext.mode ?: 'run'
    def database_command = database ? "--amrfinder_db ${database}" : ""
    def matches_command = task.ext.mode == 'report' ? "--matches ${matches}" : ""
    def partials_command = task.ext.mode == 'report' ? "--partials ${partials}" : ""
    def args = task.ext.args ?: ''
    matches = "${sampleID}_matches.tsv"
    partials = "${sampleID}_partials.tsv"
    virulence = "${sampleID}_virulence.tsv"
    output = "${sampleID}_amr.out"
    """
    abritamr $mode \\
    --contigs $contigs \\
    --species $species \\
    $database_command \\
    $matches_command \\
    $partials_command \\
    $args

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     abritamr:
      version: \$(echo \$(abritamr --version 2>&1) | sed 's/^.*abritamr // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    """
    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     abritamr:
      version: \$(echo \$(abritamr --version 2>&1) | sed 's/^.*abritamr // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
