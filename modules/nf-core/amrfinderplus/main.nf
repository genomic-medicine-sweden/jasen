include { get_species_taxon_name } from '../../../methods/get_taxon.nf'

process amrfinderplus {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(assembly)
    val species
    path database

  output:
    tuple val(sample_id), path(output), emit: tsv
    path "*versions.yml"              , emit: versions

  when:
    task.ext.when

  script:
    def args = task.ext.args ?: ''
    def database_arg = database ? "--database ${database}" : ""
    def taxon_name = get_species_taxon_name(species)
    def taxon_arg = taxon_name ? "--organism ${taxon_name}" : ""
    output = "${sample_id}_amrfinder.tsv"
    """
    amrfinder \\
    --nucleotide ${assembly} \\
    ${database_arg} \\
    ${args} \\
    ${taxon_arg} \\
    --output ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     amrfinderplus:
      version: \$(echo \$(amrfinder --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sample_id}_amrfinder.tsv"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     amrfinderplus:
      version: \$(echo \$(amrfinder --version 2>&1))
      container: ${task.container}
    END_VERSIONS
    """
}
