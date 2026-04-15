#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { clair3    } from '../modules/nf-core/clair3/main.nf'
include { freebayes } from '../modules/nf-core/freebayes/main.nf'

workflow CALL_VARIANT_CALLING {
    take:
    ch_ref_bam
    ch_ref_bai
    reference_genome
    reference_genome_faidx
    clair3_model

    main:

    ch_versions = Channel.empty()

    ch_ref_bam
        .join(ch_ref_bai)
        .set{ ch_ref_bam_bai }

    if ( params.platform == "nanopore" ) {
        clair3(ch_ref_bam_bai, reference_genome, reference_genome_faidx, clair3_model)
        ch_vcf = clair3.out.vcf
        ch_versions = ch_versions.mix(clair3.out.versions)
    } else {
        freebayes(
            ch_ref_bam_bai.map { id, bam, bai -> tuple(id, reference_genome, bam, bai) }
        )
        ch_vcf = freebayes.out.vcf
        ch_versions = ch_versions.mix(freebayes.out.versions)
    }

    emit:
    vcf         = ch_vcf        // channel: [ val(meta), path(vcf) ]
    versions    = ch_versions   // channel: [ versions.yml ]
}
