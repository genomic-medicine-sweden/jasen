include { get_species_taxon_name } from '../../../methods/get_taxon.nf'

process amrfinderplus {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(assembly)
    val species
    path database

  output:
    tuple val(sampleID), path(output), emit: output
    path "*versions.yml"             , emit: versions

  script:
    def args = task.ext.args ?: ''
    def database_command = database ? "--database ${database}" : ""
    def taxon_name = get_species_taxon_name(species)
    def taxon_command = taxon_name ? "--organism $taxon_name" : ""
    output = "${sampleID}_amrfinder.out"
    """
    amrfinder \\
    --nucleotide $assembly \\
    $database_command \\
    $args \\
    $taxon_command \\
    --output $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     amrfinderplus:
      version: \$(echo \$(amrfinder --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleID}_amrfinder.out"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     amrfinderplus:
      version: \$(echo \$(amrfinder --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """
}
