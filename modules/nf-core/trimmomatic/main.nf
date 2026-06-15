process trimmomatic {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path(output), emit: reads
    path "*versions.yml"              , emit: versions
    path "*.log"                      , emit: log

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    def is_paired = (reads instanceof List) && reads.size() == 2
    if (is_paired) {
        output = [
            "${sample_id}.paired.trim_1.fastq.gz",
            "${sample_id}.paired.trim_2.fastq.gz",
        ]
        """
        trimmomatic PE \\
            -threads ${task.cpus} \\
            ${reads[0]} ${reads[1]} \\
            ${sample_id}.paired.trim_1.fastq.gz ${sample_id}.unpaired.trim_1.fastq.gz \\
            ${sample_id}.paired.trim_2.fastq.gz ${sample_id}.unpaired.trim_2.fastq.gz \\
            ${args} \\
            2> >(tee ${sample_id}.trimmomatic.log >&2)

        cat <<END_VERSIONS > ${sample_id}_${task.process}_versions.yml
        ${task.process}:
         trimmomatic:
          version: \$(echo \$(trimmomatic -version 2>&1))
          container: ${task.container}
END_VERSIONS
        """
    } else {
        def single_input = (reads instanceof List) ? reads[0] : reads
        output = "${sample_id}.SE.trim.fastq.gz"
        """
        trimmomatic SE \\
            -threads ${task.cpus} \\
            ${single_input} \\
            ${output} \\
            ${args} \\
            2> >(tee ${sample_id}.trimmomatic.log >&2)

        cat <<END_VERSIONS > ${sample_id}_${task.process}_versions.yml
        ${task.process}:
         trimmomatic:
          version: \$(echo \$(trimmomatic -version 2>&1))
          container: ${task.container}
END_VERSIONS
        """
    }

    stub:
    def is_paired = (reads instanceof List) && reads.size() == 2
    if (is_paired) {
        output = [
            "${sample_id}.paired.trim_1.fastq.gz",
            "${sample_id}.paired.trim_2.fastq.gz",
        ]
        """
        touch ${sample_id}.paired.trim_1.fastq.gz
        touch ${sample_id}.paired.trim_2.fastq.gz
        touch ${sample_id}.trimmomatic.log

        cat <<END_VERSIONS > ${sample_id}_${task.process}_versions.yml
        ${task.process}:
         trimmomatic:
          version: \$(echo \$(trimmomatic -version 2>&1))
          container: ${task.container}
END_VERSIONS
        """
    } else {
        output = "${sample_id}.SE.trim.fastq.gz"
        """
        touch ${output}
        touch ${sample_id}.trimmomatic.log

        cat <<END_VERSIONS > ${sample_id}_${task.process}_versions.yml
        ${task.process}:
         trimmomatic:
          version: \$(echo \$(trimmomatic -version 2>&1))
          container: ${task.container}
END_VERSIONS
        """
    }
}
