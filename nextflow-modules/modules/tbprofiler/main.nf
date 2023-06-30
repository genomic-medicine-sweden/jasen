process tbprofiler {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(reads)

  output:
    tuple val(sampleName), path("results/*.json"), emit: json
    path "*versions.yml"                         , emit: versions

  when:
    task.ext.when

  script:
    def args = task.ext.args ?: ''
    def inputData = reads.size() == 2 ? "-1 ${reads[0]} -2 ${reads[1]}" : "-1 ${reads[0]}"
    output = "${sampleName}.json"
    """
    tb-profiler profile \\
      ${args} \\
      ${inputData} \\
      --threads ${task.cpus} \\
      --ram ${task.memory} \\
      --prefix ${sampleName}

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     tbprofiler:
      version: \$(echo \$(tb-profiler version 2>&1) | sed 's/^.*TBProfiler version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    """
    mkdir results
    touch results/${sampleName}.json

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     tbprofiler:
      version: \$(echo \$(tb-profiler version 2>&1) | sed 's/^.*TBProfiler version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
