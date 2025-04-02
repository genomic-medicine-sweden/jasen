#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ska_build                                  } from '../modules/nf-core/ska/main.nf'
include { sourmash                                   } from '../modules/nf-core/sourmash/main.nf'

workflow CALL_RELATEDNESS {
    take:
    ch_assembly
    ch_reads

    main:

    ch_versions = Channel.empty()

    // RELATEDNESS
    sourmash(ch_assembly)

    ska_build(ch_reads)

    ska_build.out.skf
        .join(sourmash.out.signature)
        .set{ ch_combined_output }

    ch_versions = ch_versions.mix(ska_build.out.versions)
    ch_versions = ch_versions.mix(sourmash.out.versions)

    emit:
    combined_output = ch_combined_output        // channel: [ val(meta), path(skf), path(signature) ]
    ska             = ska_build.out.skf         // channel: [ val(meta), path(skf) ]
    sourmash        = sourmash.out.signature    // channel: [ val(meta), path(signature) ]
    versions        = ch_versions               // channel: [ versions.yml ]
}
