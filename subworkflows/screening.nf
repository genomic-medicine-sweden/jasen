include { amrfinderplus     } from '../nextflow-modules/modules/amrfinderplus/main'
include { resfinder         } from '../nextflow-modules/modules/resfinder/main'
include { virulencefinder   } from '../nextflow-modules/modules/virulencefinder/main'

workflow CALL_SCREENING {
    take:
        amrfinderDb
        pointfinderDb
        resfinderDb
        virulencefinderDb
        ch_assembly         // channel: [ val(meta), val(fasta) ]
        ch_reads            // channel: [ val(meta), val(reads) ]

    main:
        ch_versions = Channel.empty()

        // antimicrobial detection (amrfinderplus & abritamr)
        amrfinderplus(ch_assembly, amrfinderDb)
        //abritamr(amrfinderplus.out.output)

        // perform resistance prediction
        resfinder(ch_reads, params.species, resfinderDb, pointfinderDb)
        virulencefinder(ch_reads, params.useVirulenceDbs, virulencefinderDb)

        ch_versions = ch_versions.mix(amrfinderplus.out.versions)
        ch_versions = ch_versions.mix(resfinder.out.versions)
        ch_versions = ch_versions.mix(virulencefinder.out.versions)

    emit:
        amrfinderplus       = amrfinderplus.out.output  // channel: [ val(meta), path(tsv)]
        resfinderJson       = resfinder.out.json        // channel: [ val(meta), path(json)]
        resfinderMeta       = resfinder.out.meta        // channel: [ val(meta), path(meta)]
        virulencefinderJson = virulencefinder.out.json  // channel: [ val(meta), path(json)]
        virulencefinderMeta = virulencefinder.out.meta  // channel: [ val(meta), path(meta)]
        versions            = ch_versions               // channel: [ versions.yml ]
}