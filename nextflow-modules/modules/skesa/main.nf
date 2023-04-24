process skesa {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(reads), val(platform) 

  output:
    tuple val(sampleName), path(output), emit: fasta
    path "*versions.yml"                    , emit: versions

  when:
    task.ext.when && platform == "illumina"

  script:
    def args = task.ext.args ?: ''
    def inputData = reads.size() == 2 ? "${reads[0]},${reads[1]}" : "${reads[0]}"
    output = "${sampleName}.fasta"
    """
    skesa --reads ${inputData} ${args} > ${output}

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     skesa:
      version: \$(echo \$(skesa --version 2>&1) | sed 's/^.*SKESA // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleName}.fasta"
    """
    touch $output

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     skesa:
      version: \$(echo \$(skesa --version 2>&1) | sed 's/^.*SKESA // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
