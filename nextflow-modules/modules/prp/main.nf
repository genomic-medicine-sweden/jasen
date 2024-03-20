process create_analysis_result {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(quast), path(postalignqc), path(mlst), path(cgmlst), path(amr), path(resistance), path(resfinderMeta), path(serotype), path(serotypefinderMeta), path(virulence), path(virulencefinderMeta), path(runInfo), path(dellyVcf), path(mykrobe), path(tbprofiler), path(bracken)

  output:
    tuple val(sampleName), path(output), emit: json
    path "*versions.yml"               , emit: versions

  script:
    output = "${sampleName}_result.json"
    amrfinderArgs = amr ? "--amrfinder ${amr}" : ""
    brackenArgs = bracken ? "--kraken ${bracken}" : ""
    cgmlstArgs = cgmlst ? "--cgmlst ${cgmlst}" : ""
    dellyVcfArgs = dellyVcf ? "--sv-vcf ${dellyVcf}" : ""
    mlstArgs = mlst ? "--mlst ${mlst}" : ""
    mykrobeArgs = mykrobe ? "--mykrobe ${mykrobe}" : ""
    postalignqcArgs = postalignqc ? "--quality ${postalignqc}" : "" 
    quastArgs = quast ? "--quast ${quast}" : ""
    resfinderArgs = resistance ? "--resfinder ${resistance}" : ""
    resfinderArgs = resfinderMeta ? "${resfinderArgs} --process-metadata ${resfinderMeta}" : resfinderArgs
    runInfoArgs = runInfo ? "--run-metadata ${runInfo}" : ""
    serotypeArgs = serotype ? "--serotypefinder ${serotype}" : ""
    serotypeArgs = serotypefinderMeta ? "${serotypeArgs} --process-metadata ${serotypefinderMeta}" : serotypeArgs
    tbprofilerArgs = tbprofiler ? "--tbprofiler ${tbprofiler}" : ""
    virulenceArgs = virulence ? "--virulencefinder ${virulence}" : ""
    virulenceArgs = virulencefinderMeta ? "${virulenceArgs} --process-metadata ${virulencefinderMeta}" : virulenceArgs
    """
    prp create-bonsai-input \\
      --sample-id ${sampleName} \\
      ${amrfinderArgs} \\
      ${brackenArgs} \\
      ${cgmlstArgs} \\
      ${dellyVcfArgs} \\
      ${mlstArgs} \\
      ${mykrobeArgs} \\
      ${postalignqcArgs} \\
      ${quastArgs} \\
      ${resfinderArgs} \\
      ${runInfoArgs} \\
      ${serotypeArgs} \\
      ${tbprofilerArgs} \\
      ${virulenceArgs} \\
      --output ${output}

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     prp:
      version: \$(echo \$(prp --version 2>&1) | sed 's/prp, version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """

  stub:
    output = "${sampleName}_result.json"
    """
    touch $output

    cat <<-END_VERSIONS > ${sampleName}_${task.process}_versions.yml
    ${task.process}:
     prp:
      version: \$(echo \$(prp --version 2>&1) | sed 's/prp, version // ; s/ .*//')
      container: ${task.container}
    END_VERSIONS
    """
}

process create_cdm_input {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), val(quast), val(postalignqc), val(cgmlst)

  output:
    tuple val(sampleName), path(output), emit: json

  script:
    output = "${sampleName}_qc_result.json"
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
    output = "${sampleName}_qc_result.json"
    """
    touch $output
    """
}

process post_align_qc {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(bam)
    path reference
    path bed

  output:
    tuple val(sampleName), path(output), emit: qc

  script:
    output = "${sampleName}_bwa.qc"
    """
    prp create-qc-result --bam ${bam} --reference ${reference} --bed ${bed} --sample-id ${sampleName} --cpus ${task.cpus} --output ${output}
    """

  stub:
    output = "${sampleName}_bwa.qc"
    """
    touch $output
    """
}

process annotate_delly {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(vcf)
    path bed

  output:
    tuple val(sampleName), path(output), emit: vcf

  script:
    output = "${sampleName}_annotated_delly.vcf"
    """
    prp annotate-delly --vcf ${vcf} --bed ${bed} --output ${output}
    """

  stub:
    output = "${sampleName}_annotated_delly.vcf"
    """
    touch $output
    """
}
