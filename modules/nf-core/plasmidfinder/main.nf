process plasmidfinder {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(assembly)
    path plasmidfinder_db

    output:
    tuple val(sample_id), path(output)     , emit: json
    tuple val(sample_id), path(meta_output), emit: meta
    tuple val(sample_id), path("${sample_id}_plasmidfinder_hit_in_genome_seq.fsa"), emit: genome_hits
    tuple val(sample_id), path("${sample_id}_plasmidfinder_plasmid_seqs.fsa")    , emit: plasmid_seqs
    tuple val(sample_id), path("${sample_id}_plasmidfinder_results.tsv"), optional: true, emit: tsv
    tuple val(sample_id), path("${sample_id}_plasmidfinder_results.txt"), optional: true, emit: txt
    path "*versions.yml"                   , emit: versions

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    output = "${sample_id}_plasmidfinder.json"
    meta_output = "${sample_id}_plasmidfinder_meta.json"
    """
    # Get db version
    if [ -f "${plasmidfinder_db}/VERSION" ]; then
        DB_VERSION=\$(tr -d '\r\n' < ${plasmidfinder_db}/VERSION)
    else
        DB_VERSION=\$(cd ${plasmidfinder_db} && (git rev-parse --short HEAD 2>/dev/null || echo "unknown"))
    fi
    JSON_FMT='{"name": "%s", "version": "%s", "type": "%s"}'
    printf "\$JSON_FMT" "plasmidfinder" "\$DB_VERSION" "database" > ${meta_output}

    # Run plasmidfinder
    plasmidfinder.py \\
        ${args} \\
        -i ${assembly} \\
        -o ./ \\
        -p ${plasmidfinder_db} \\
        -x

    # Rename hard-coded outputs to sample-specific names
    mv data.json ${output}
    mv Hit_in_genome_seq.fsa ${sample_id}_plasmidfinder_hit_in_genome_seq.fsa
    mv Plasmid_seqs.fsa ${sample_id}_plasmidfinder_plasmid_seqs.fsa
    [ -f results.txt ] && mv results.txt ${sample_id}_plasmidfinder_results.txt || true
    [ -f results_tab.tsv ] && mv results_tab.tsv ${sample_id}_plasmidfinder_results.tsv || true

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     plasmidfinder:
      version: \$(echo \$(plasmidfinder.py --version 2>&1) | sed 's/^.*plasmidfinder.py //; s/ .*\$//')
      container: ${task.container}
     plasmidfinder_db:
      version: \$(echo \$DB_VERSION)
      container: ${task.container}
    END_VERSIONS
    """

    stub:
    output = "${sample_id}_plasmidfinder.json"
    meta_output = "${sample_id}_plasmidfinder_meta.json"
    """
    touch ${output}
    touch ${meta_output}
    touch ${sample_id}_plasmidfinder_hit_in_genome_seq.fsa
    touch ${sample_id}_plasmidfinder_plasmid_seqs.fsa

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     plasmidfinder:
      version: stub
      container: ${task.container}
    END_VERSIONS
    """
}
