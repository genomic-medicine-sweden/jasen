process format_jasen {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(yaml)

    output:
    tuple val(sample_id), path(output), emit: json
    path "*versions.yml"              , emit: versions

    script:
    output = "${sample_id}_result.json"
    """
    prp format-jasen \\
      --sample ${yaml} \\
      --output ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     prp:
      version: \$(echo \$(prp --version 2>&1) | sed 's/prp, version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

    stub:
    output = "${sample_id}_result.json"
    """
    touch ${output}

    cat <<-END_VERSIONS > ${sample_id}_${task.process}_versions.yml
    ${task.process}:
     prp:
      version: \$(echo \$(prp --version 2>&1) | sed 's/prp, version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}

process format_cdm {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(yaml)

    output:
    tuple val(sample_id), path(output), emit: json

    script:
    output = "${sample_id}_qc_result.json"
    """
    prp format-cdm \\
      --sample ${yaml} \\
      --output ${output}
    """

    stub:
    output = "${sample_id}_qc_result.json"
    """
    touch ${output}
    """
}

process post_align_qc {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(bam)
    path reference
    path bed

    output:
    tuple val(sample_id), path(output), emit: json

    when:
    task.ext.when

    script:
    output = "${sample_id}_qc.json"
    """
    prp create-qc-result --bam ${bam} --reference ${reference} --bed ${bed} --sample-id ${sample_id} --cpus ${task.cpus} --output ${output}
    """

    stub:
    output = "${sample_id}_qc.json"
    """
    touch ${output}
    """
}

process annotate_delly {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(vcf)
    path bed
    path bedIdx

    output:
    tuple val(sample_id), path(output), emit: vcf

    when:
    task.ext.when

    script:
    output = "${sample_id}_annotated_delly.vcf"
    """
    prp annotate-delly --vcf ${vcf} --bed ${bed} --output ${output}
    """

    stub:
    output = "${sample_id}_annotated_delly.vcf"
    """
    touch ${output}
    """
}

process add_igv_track {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(bonsaiInput)
    val annotation
    val trackName

    output:
    tuple val(sample_id), path(output), emit: json

    when:
    task.ext.when

    script:
    output = "${sample_id}_result.json"
    """
    prp add-igv-annotation-track --track-name ${trackName} --annotation-file ${annotation} --bonsai-input-file ${bonsaiInput} --output ${output}
    """

    stub:
    output = "${sample_id}_result.json"
    """
    touch ${output}
    """
}
