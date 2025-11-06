include { check_taxon } from '../../../methods/check_taxon.nf'

process kleborate {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("*.txt")    , emit: txt
    path "*versions.yml"                   , emit: versions

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Run kleborate
    kleborate                \\
    ${args}                  \\
    --outdir results         \\
    --assemblies ${assembly}

    # Move results to cwd
    mv results/*.txt .

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     kleborate:
      version: \$(echo \$(kleborate --version 3>&1 | sed "s/.*v//") )
      container: ${task.container}
    END_VERSIONS
    """

    stub:
    """
    touch results.txt

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     kleborate:
      version: \$(echo \$(kleborate --version 2>&1 | sed "s/.*v//") )
      container: ${task.container}
    END_VERSIONS
    """
}
