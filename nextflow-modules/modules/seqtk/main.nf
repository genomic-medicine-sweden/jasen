process seqtk_sample {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(reads), val(sample_size)

  output:
    tuple val(sampleID), path("*.fastq.gz"), emit: reads
    path "*versions.yml"                   , emit: versions

  when:
    task.ext.when

  script:
    def args   = task.ext.args ?: ''
    if (!(args ==~ /.*-s[0-9]+.*/)) {
        args += " -s100"
    }
    if ( !sample_size ) {
        error "SEQTK/SAMPLE must have a sample_size value included"
    }
    """
<<<<<<< HEAD
    printf "%s\\n" $reads | while read f;
    do
=======
    printf "%s\\n" $reads | while read f; do
>>>>>>> c5624dc39bf9062bd570a19a80548ebc8ee0eb3a
        output_name=\$(basename \$f | sed "s/.fastq/_seqtk.fastq/")
        seqtk \\
            sample \\
            $args \\
            \$f \\
            $sample_size \\
            | gzip --no-name > \${output_name}
    done

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     seqtk:
      version: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleID}_seqtk.fastq.gz"
    """
    touch $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     seqtk:
      version: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
      container: ${task.container}
    END_VERSIONS
    """
}
