process create_analysis_result {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(quast), path(postalignqc), path(mlst), path(cgmlst), path(amr), path(resistance), path(resfinderMeta), path(serotype), path(serotypefinderMeta), path(virulence), path(virulencefinderMeta), path(shigapass), path(bam), path(bai), path(runInfo), path(vcf), path(mykrobe), path(tbprofiler), path(bracken)
    path referenceGenome
    path referenceGenomeIdx
    path referenceGenomeGff

  output:
    tuple val(sampleID), path(output), emit: json
    path "*versions.yml"             , emit: versions

  script:
    output = "${sampleID}_result.json"
    amrfinderArgs = amr ? "--amrfinder ${amr}" : ""
    brackenArgs = bracken ? "--kraken ${bracken}" : ""
    bamArgs = bam ? "--bam ${params.outdir}/${params.speciesDir}/${params.bamDir}/${bam}" : ""
    cgmlstArgs = cgmlst ? "--cgmlst ${cgmlst}" : ""
    vcfArgs = vcf ? "--vcf ${params.outdir}/${params.speciesDir}/${params.vcfDir}/${vcf}" : ""
    mlstArgs = mlst ? "--mlst ${mlst}" : ""
    mykrobeArgs = mykrobe ? "--mykrobe ${mykrobe}" : ""
    postalignqcArgs = postalignqc ? "--quality ${postalignqc}" : "" 
    quastArgs = quast ? "--quast ${quast}" : ""
    referenceGenomeArgs = referenceGenome ? "--reference-genome-fasta ${referenceGenome}" : ""
    referenceGffArgs = referenceGenomeGff ? "--reference-genome-gff ${referenceGenomeGff}" : ""
    resfinderArgs = resistance ? "--resfinder ${resistance}" : ""
    resfinderArgs = resfinderMeta ? "${resfinderArgs} --process-metadata ${resfinderMeta}" : resfinderArgs
    runInfoArgs = runInfo ? "--run-metadata ${runInfo}" : ""
    serotypeArgs = serotype ? "--serotypefinder ${serotype}" : ""
    serotypeArgs = serotypefinderMeta ? "${serotypeArgs} --process-metadata ${serotypefinderMeta}" : serotypeArgs
    shigapassArgs = shigapass ? "--shigapass ${shigapass}" : ""
    symlinkDirArgs = params.symlinkDir ? "--symlink-dir ${params.symlinkDir}" : ""
    tbprofilerArgs = tbprofiler ? "--tbprofiler ${tbprofiler}" : ""
    virulenceArgs = virulence ? "--virulencefinder ${virulence}" : ""
    virulenceArgs = virulencefinderMeta ? "${virulenceArgs} --process-metadata ${virulencefinderMeta}" : virulenceArgs
    """
    prp create-bonsai-input \\
      --sample-id ${sampleID} \\
      ${amrfinderArgs} \\
      ${bamArgs} \\
      ${brackenArgs} \\
      ${cgmlstArgs} \\
      ${vcfArgs} \\
      ${mlstArgs} \\
      ${mykrobeArgs} \\
      ${postalignqcArgs} \\
      ${quastArgs} \\
      ${referenceGenomeArgs} \\
      ${referenceGffArgs} \\
      ${resfinderArgs} \\
      ${runInfoArgs} \\
      ${serotypeArgs} \\
      ${shigapassArgs} \\
      ${symlinkDirArgs} \\
      ${tbprofilerArgs} \\
      ${virulenceArgs} \\
      --output ${output}

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     prp:
      version: \$(echo \$(prp --version 2>&1) | sed 's/prp, version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleID}_result.json"
    """
    touch $output

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     prp:
      version: \$(echo \$(prp --version 2>&1) | sed 's/prp, version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}

process create_cdm_input {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), val(quast), val(postalignqc), val(cgmlst)

  output:
    tuple val(sampleID), path(output), emit: json

  script:
    output = "${sampleID}_qc_result.json"
    cgmlstArgs = cgmlst ? "--cgmlst ${cgmlst}" : ""
    postalignqcArgs = postalignqc ? "--quality ${postalignqc}" : "" 
    quastArgs = quast ? "--quast ${quast}" : ""
    """
    prp create-cdm-input \\
      ${cgmlstArgs} \\
      ${postalignqcArgs} \\
      ${quastArgs} \\
      --output ${output}
    """

  stub:
    output = "${sampleID}_qc_result.json"
    """
    touch $output
    """
}

process post_align_qc {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(bam)
    path reference
    path bed

  output:
    tuple val(sampleID), path(output), emit: qc

  script:
    output = "${sampleID}_qc.json"
    """
    prp create-qc-result --bam ${bam} --reference ${reference} --bed ${bed} --sample-id ${sampleID} --cpus ${task.cpus} --output ${output}
    """

  stub:
    output = "${sampleID}_qc.json"
    """
    touch $output
    """
}

process annotate_delly {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(vcf)
    path bed
    path bedIdx

  output:
    tuple val(sampleID), path(output), emit: vcf

  script:
    output = "${sampleID}_annotated_delly.vcf"
    """
    prp annotate-delly --vcf ${vcf} --bed ${bed} --output ${output}
    """

  stub:
    output = "${sampleID}_annotated_delly.vcf"
    """
    touch $output
    """
}

process add_igv_track {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(bonsaiInput)
    val annotation
    val trackName

  output:
    tuple val(sampleID), path(output), emit: json

  script:
    output = "${sampleID}_result.json"
    """
    prp add-igv-annotation-track --track-name ${trackName} --annotation-file ${annotation} --bonsai-input-file ${bonsaiInput} --output ${output}
    """

  stub:
    output = "${sampleID}_result.json"
    """
    touch $output
    """
}
