process kraken_batch {
    tag "${workflow.runName}"
    scratch params.scratch

    input:
    path batch_input
    path database

    output:
    path("*_kraken.out"),    emit: outputs
    path("*_kraken.report"), emit: reports
    path("*versions.yml"),   emit: versions

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    SHM_DB="/dev/shm/$(basename ${database})_\${BASHPID}"

    cleanup() {
        rm -rf "\${SHM_DB}"
    }
    trap cleanup EXIT SIGINT SIGTERM

    if [[ -d "\${SHM_DB}" ]]; then
        echo "Error: \${SHM_DB} already exists. Clean up first."
        trap - EXIT
        exit 1
    fi
    mkdir -p "\${SHM_DB}"

    for file in ${database}/*; do
        if ! dd if="\${file}" of="\${SHM_DB}/\${file##*/}" \\
             bs=4M conv=fdatasync iflag=nocache oflag=nocache \\
             status=progress; then
            echo "ERROR: Failed to copy \${file##*/}"
            exit 1
        fi
    done

    while IFS=\$'\\t' read -r sample_id fastq1 fastq2; do
        reads_arg="\${fastq1}\${fastq2:+ \${fastq2}}"
        kraken2 \\
        ${args} \\
        --memory-mapping \\
        --threads ${task.cpus} \\
        --db \${SHM_DB} \\
        --output \${sample_id}_kraken.out \\
        --report \${sample_id}_kraken.report \\
        \${reads_arg}
    done < ${batch_input}

    cat <<-END_VERSIONS > kraken_batch_versions.yml
    ${task.process}:
     kraken2:
      version: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

    stub:
    """
    while IFS=\$'\\t' read -r sample_id fastq1 fastq2; do
        touch \${sample_id}_kraken.out
        touch \${sample_id}_kraken.report
    done < ${batch_input}

    cat <<-END_VERSIONS > kraken_batch_versions.yml
    ${task.process}:
     kraken2:
      version: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}

process kraken {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads)
    path database

    output:
    tuple val(sample_id), path(output), emit: output
    tuple val(sample_id), path(report), emit: report
    path "*versions.yml"              , emit: versions

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_reads_arg = reads.size() == 2 ? "${reads[0]} ${reads[1]}" : "${reads[0]}"
    output = "${sample_id}_kraken.out"
    report = "${sample_id}_kraken.report"
    """
    kraken2 \\
    ${args} \\
    --threads ${task.cpus} \\
    --db ${database} \\
    --output ${output} \\
    --report ${report} \\
    ${input_reads_arg}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     kraken2:
      version: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

    stub:
    output = "${sample_id}_kraken.out"
    report = "${sample_id}_kraken.report"
    """
    touch ${output}
    touch ${report}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     kraken2:
      version: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
