#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { freebayes } from '../modules/nf-core/freebayes/main.nf'

workflow CALL_VARIANT_CALLING {
    take:
    ch_ref_bam
    ch_ref_bai
    reference_genome

    main:

    ch_versions = Channel.empty()

    // VARIANT CALLING against reference genome (for reporting / postprocessing)
    ch_ref_bam
        .join(ch_ref_bai)
        .set{ ch_ref_bam_bai }

    freebayes(
        ch_ref_bam_bai.map { id, bam, bai -> tuple(id, reference_genome, bam, bai) }
    )
    ch_vcf = freebayes.out.vcf

    ch_versions = ch_versions.mix(freebayes.out.versions)

    emit:
    vcf         = ch_vcf                // channel: [ val(meta), path(vcf) ]
    versions    = ch_versions           // channel: [ versions.yml ]
}
