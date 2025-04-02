process seqtk_sample {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads)
    val sample_size

    output:
    tuple val(sample_id), path("*.fastq.gz"), emit: reads
    path "*versions.yml"                    , emit: versions

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    if (!(args ==~ /.*-s[0-9]+.*/)) {
        args += " -s100"
    }
    if ( !sample_size ) {
        error "SEQTK/SAMPLE must have a sample_size value included"
    }
    """
    printf "%s\\n" ${reads} | while read f; do
        output_name=\$(basename \$f | sed "s/.fastq/_seqtk.fastq/")
        seqtk \\
            sample \\
            ${args} \\
            \$f \\
            ${sample_size} \\
            | gzip --no-name > \${output_name}
    done

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     seqtk:
      version: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
      container: ${task.container}
    END_VERSIONS
    """

    stub:
    output = "${sample_id}_seqtk.fastq.gz"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     seqtk:
      version: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
      container: ${task.container}
    END_VERSIONS
    """
}
