process tbprofiler {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(reads)

  output:
    tuple val(sampleID), path(vcfOutput), emit: vcf
    tuple val(sampleID), path(output)   , emit: json
    tuple val(sampleID), path(bamOutput), emit: bam
    tuple val(sampleID), path(baiOutput), emit: bai
    path "*versions.yml"                , emit: versions

  when:
    workflow.profile == "mtuberculosis"

  script:
    def args = task.ext.args ?: ''
    def inputData = reads.size() == 2 ? "-1 ${reads[0]} -2 ${reads[1]}" : "-1 ${reads[0]}"
    output = "${sampleID}_tbprofiler.json"
    vcfOutput = "${sampleID}_tbprofiler.vcf.gz"
    bamOutput = "${sampleID}_tbprofiler.bam"
    baiOutput = "${bamOutput}.bai"
    """
    tb-profiler profile \\
      ${args} \\
      ${inputData} \\
      --threads ${task.cpus} \\
      --prefix ${sampleID}

    cp results/${sampleID}.results.json $output
    cp vcf/${sampleID}.targets.vcf.gz $vcfOutput
    cp bam/${sampleID}.bam $bamOutput
    cp bam/${sampleID}.bam.bai $baiOutput

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     tbprofiler:
      version: \$(echo \$(tb-profiler version 2>&1) | sed 's/^.*TBProfiler version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleID}_tbprofiler.json"
    vcfOutput = "${sampleID}_tbprofiler.vcf.gz"
    bamOutput = "${sampleID}_tbprofiler.bam"
    baiOutput = "${bamOutput}.bai"
    """
    mkdir results
    touch $output
    touch $vcfOutput
    touch $bamOutput
    touch $baiOutput

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     tbprofiler:
      version: \$(echo \$(tb-profiler version 2>&1) | sed 's/^.*TBProfiler version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
