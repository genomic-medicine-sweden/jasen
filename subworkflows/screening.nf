#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { amrfinderplus     } from '../modules/nf-core/amrfinderplus/main.nf'
include { resfinder         } from '../modules/nf-core/resfinder/main.nf'
include { virulencefinder   } from '../modules/nf-core/virulencefinder/main.nf'

workflow CALL_SCREENING {
    take:
    amrfinder_db
    pointfinder_db
    resfinder_db
    virulencefinder_db
    ch_assembly
    ch_reads
    ch_reads_w_meta

    main:

    ch_versions = Channel.empty()

    // SCREENING
    // antimicrobial detection (amrfinderplus)
    amrfinderplus(ch_assembly, params.species, amrfinder_db)

    // resistance & virulence prediction
    resfinder(ch_reads_w_meta, params.species, resfinder_db, pointfinder_db)
    virulencefinder(ch_reads, params.use_virulence_dbs, virulencefinder_db)

    amrfinderplus.out.tsv
        .join(resfinder.out.json)
        .join(resfinder.out.meta)
        .join(virulencefinder.out.json)
        .join(virulencefinder.out.meta)
        .set{ ch_combined_output }

    ch_versions = ch_versions.mix(amrfinderplus.out.versions)
    ch_versions = ch_versions.mix(resfinder.out.versions)
    ch_versions = ch_versions.mix(virulencefinder.out.versions)

    emit:
    amrfinderplus           = amrfinderplus.out.tsv     // channel: [ val(meta), path(tsv) ]
    combined_output         = ch_combined_output        // channel: [ val(meta), path(tsv) ]
    resfinder_json          = resfinder.out.json        // channel: [ val(meta), path(json) ]
    resfinder_meta          = resfinder.out.meta        // channel: [ val(meta), path(meta) ]
    virulencefinder_json    = virulencefinder.out.json  // channel: [ val(meta), path(json) ]
    virulencefinder_meta    = virulencefinder.out.meta  // channel: [ val(meta), path(meta) ]
    versions                = ch_versions               // channel: [ versions.yml ]
}
