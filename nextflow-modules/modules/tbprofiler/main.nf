process tbprofiler {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path(vcf_output) , emit: vcf
    tuple val(sample_id), path(output)     , emit: json
    tuple val(sample_id), path(bam_output) , emit: bam
    tuple val(sample_id), path(bai_output) , emit: bai
    path "*versions.yml"                   , emit: versions

  when:
    workflow.profile == "mycobacterium_tuberculosis"

  script:
    def args = task.ext.args ?: ''
    def input_reads_arg = reads.size() == 2 ? "-1 ${reads[0]} -2 ${reads[1]}" : "-1 ${reads[0]}"
    output = "${sample_id}_tbprofiler.json"
    vcf_output = "${sample_id}_tbprofiler.vcf.gz"
    bam_output = "${sample_id}_tbprofiler.bam"
    bai_output = "${bam_output}.bai"
    """
    tb-profiler profile \\
      ${args} \\
      ${input_reads_arg} \\
      --threads ${task.cpus} \\
      --prefix ${sample_id}

    cp results/${sample_id}.results.json ${output}
    cp vcf/${sample_id}.targets.vcf.gz ${vcf_output}
    cp bam/${sample_id}.bam ${bam_output}
    cp bam/${sample_id}.bam.bai ${bai_output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     tbprofiler:
      version: \$(echo \$(tb-profiler version 2>&1) | sed 's/^.*TBProfiler version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sample_id}_tbprofiler.json"
    vcf_output = "${sample_id}_tbprofiler.vcf.gz"
    bam_output = "${sample_id}_tbprofiler.bam"
    bai_output = "${bam_output}.bai"
    """
    mkdir results
    touch ${output}
    touch ${vcf_output}
    touch ${bam_output}
    touch ${bai_output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     tbprofiler:
      version: \$(echo \$(tb-profiler version 2>&1) | sed 's/^.*TBProfiler version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}
