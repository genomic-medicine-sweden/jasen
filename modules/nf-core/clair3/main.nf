process clair3 {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path reference_fai
    val model

    output:
    tuple val(sample_id), path("${sample_id}_clair3.vcf.gz"), emit: vcf
    path "*versions.yml"                                     , emit: versions

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    def platform = task.ext.platform ?: 'ont'
    def model_path = model.toString().startsWith('/') ? model : "/opt/models/${model}"
    output = "${sample_id}_clair3.vcf.gz"
    """
    run_clair3.sh \\
        --bam_fn=${bam} \\
        --ref_fn=${reference} \\
        --threads=${task.cpus} \\
        --platform=${platform} \\
        --model_path=${model_path} \\
        --output=clair3_output \\
        ${args}

    cp clair3_output/merge_output.vcf.gz ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     clair3:
      version: \$(run_clair3.sh --version 2>&1 | sed 's/Clair3 v//')
      container: ${task.container}
    END_VERSIONS
    """

    stub:
    output = "${sample_id}_clair3.vcf.gz"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     clair3:
      version: \$(run_clair3.sh --version 2>&1 | sed 's/Clair3 v//')
      container: ${task.container}
    END_VERSIONS
    """
}
