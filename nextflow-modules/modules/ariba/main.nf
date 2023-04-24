process ariba_prepareref {
  tag "${fasta.simpleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    path fasta
    path metadata

  output:
    path(outputDir)

  script:
    def args = task.ext.args ?: ''
    def outputDir = "ariba_reference"
    def metadata = metadata ? "--metadata ${metadata}" : "--all_coding"
    """
    ariba prepareref \\
    --threads ${task.cpus} \\
    ${args} \\
    ${metadata} \\
    --fasta ${fasta} \\
    ${outputDir}
    """
}

process ariba_run {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(reads)
    path referenceDir

  output:
    tuple val(sampleName), path("${sampleName}_ariba_report.tsv")
    path "*versions.yml"

  script:
    def args = task.ext.args ?: ''
    outputName = "ariba_output"
    """
    ariba run ${args} --force --threads ${task.cpus} ${referenceDir} ${reads.join(' ')} ${outputName}
    cp ${outputName}/report.tsv ${sampleName}_ariba_report.tsv

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     ariba:
      version: \$(echo \$(ariba version -v 2>&1) | sed 's/.*ARIBA version: // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    """
    touch ${sampleName}_ariba_report.tsv

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     ariba:
      version: \$(echo \$(ariba version -v 2>&1) | sed 's/.*ARIBA version: // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}

process ariba_summary {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(report)
    path reads
    path "*versions.yml"

  output:
    tuple val(sampleName), path("${outputPrefix}.csv")

  script:
    def args = task.ext.args ?: ''
    outputPrefix = params.prefix ?: report.simpleName.replaceFirst('report', 'summary')
    """
    ariba summary ${args} ${outputPrefix} ${report}

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     ariba:
      version: \$(echo \$(ariba version -v 2>&1) | sed 's/.*ARIBA version: // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    outputPrefix = params.prefix ?: report.simpleName.replaceFirst('report', 'summary')
    """
    touch ${outputPrefix}.csv

    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
     ariba:
      version: \$(echo \$(ariba version -v 2>&1) | sed 's/.*ARIBA version: // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}

process ariba_summary_to_json {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(report), path(summary)
    path reference
    path reads

  output:
    tuple val(sampleName), path(output), emit: output

  script:
    output = "${summary.simpleName}_export.json"
    """
    ariba2json.pl ${reference} ${summary} ${report} > ${output}
    """

  stub:
    output = "${summary.simpleName}_export.json"
    """
    touch $output
    """
}