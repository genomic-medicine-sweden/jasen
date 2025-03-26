process create_analysis_result {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), path(quast), path(postalignqc), path(mlst), path(cgmlst), path(amr), path(resistance), path(resfinderMeta), path(serotype), path(serotypefinderMeta), path(virulence), path(virulencefinderMeta), path(shigapass), path(emmtyper), path(bam), path(bai), path(runInfo), path(vcf), path(mykrobe), path(tbprofiler), path(bracken)
    path reference_genome
    path reference_genome_idx
    path reference_genome_gff

  output:
    tuple val(sample_id), path(output), emit: json
    path "*versions.yml"              , emit: versions

  script:
    output = "${sample_id}_result.json"
    amrfinder_arg = amr ? "--amrfinder ${amr}" : ""
    bracken_arg = bracken ? "--kraken ${bracken}" : ""
    bam_arg = bam ? "--bam ${params.outdir}/${params.speciesDir}/${params.bamDir}/${bam}" : ""
    cgmlst_arg = cgmlst ? "--cgmlst ${cgmlst}" : ""
    emmtyper_arg = emmtyper ? "--emmtyper ${emmtyper}" : ""
    mlst_arg = mlst ? "--mlst ${mlst}" : ""
    mykrobe_arg = mykrobe ? "--mykrobe ${mykrobe}" : ""
    postalignqc_arg = postalignqc ? "--quality ${postalignqc}" : "" 
    quast_arg = quast ? "--quast ${quast}" : ""
    reference_genome_arg = reference_genome ? "--reference-genome-fasta ${reference_genome}" : ""
    reference_gff_arg = reference_genome_gff ? "--reference-genome-gff ${reference_genome_gff}" : ""
    resfinder_arg = resistance ? "--resfinder ${resistance}" : ""
    resfinder_arg = resfinderMeta ? "${resfinder_arg} --process-metadata ${resfinderMeta}" : resfinder_arg
    runInfo_arg = runInfo ? "--run-metadata ${runInfo}" : ""
    serotype_arg = serotype ? "--serotypefinder ${serotype}" : ""
    serotype_arg = serotypefinderMeta ? "${serotype_arg} --process-metadata ${serotypefinderMeta}" : serotype_arg
    shigapass_arg = shigapass ? "--shigapass ${shigapass}" : ""
    symlink_dir_arg = params.symlinkDir ? "--symlink-dir ${params.symlinkDir}" : ""
    tbprofiler_arg = tbprofiler ? "--tbprofiler ${tbprofiler}" : ""
    vcf_arg = vcf ? "--vcf ${params.outdir}/${params.speciesDir}/${params.vcfDir}/${vcf}" : ""
    virulence_arg = virulence ? "--virulencefinder ${virulence}" : ""
    virulence_arg = virulencefinderMeta ? "${virulence_arg} --process-metadata ${virulencefinderMeta}" : virulence_arg
    """
    prp create-bonsai-input \\
      --sample-id ${sample_id} \\
      ${amrfinder_arg} \\
      ${bam_arg} \\
      ${bracken_arg} \\
      ${cgmlst_arg} \\
      ${emmtyper_arg} \\
      ${mlst_arg} \\
      ${mykrobe_arg} \\
      ${postalignqc_arg} \\
      ${quast_arg} \\
      ${reference_genome_arg} \\
      ${reference_gff_arg} \\
      ${resfinder_arg} \\
      ${runInfo_arg} \\
      ${serotype_arg} \\
      ${shigapass_arg} \\
      ${symlink_dir_arg} \\
      ${tbprofiler_arg} \\
      ${vcf_arg} \\
      ${virulence_arg} \\
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

process create_cdm_input {
  tag "${sample_id}"
  scratch params.scratch

  input:
    tuple val(sample_id), val(quast), val(postalignqc), val(cgmlst)

  output:
    tuple val(sample_id), path(output), emit: json

  script:
    output = "${sample_id}_qc_result.json"
    cgmlst_arg = cgmlst ? "--cgmlst ${cgmlst}" : ""
    postalignqc_arg = postalignqc ? "--quality ${postalignqc}" : "" 
    quast_arg = quast ? "--quast ${quast}" : ""
    """
    prp create-cdm-input \\
      ${cgmlst_arg} \\
      ${postalignqc_arg} \\
      ${quast_arg} \\
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
    tuple val(sample_id), path(output), emit: qc

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
