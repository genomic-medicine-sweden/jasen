#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { flye                          } from '../modules/nf-core/flye/main.nf'
include { medaka                        } from '../modules/nf-core/medaka/main.nf'
include { skesa                         } from '../modules/nf-core/skesa/main.nf'
include { spades as spades_illumina     } from '../modules/nf-core/spades/main.nf'
include { spades as spades_iontorrent   } from '../modules/nf-core/spades/main.nf'

workflow CALL_ASSEMBLY {
    take:
    ch_reads

    main:

    ch_versions = Channel.empty()

    // ASSEMBLY
    skesa(ch_reads)
    spades_illumina(ch_reads)
    spades_iontorrent(ch_reads)
    flye(ch_reads)
    medaka(ch_reads, flye.out.fasta)

    Channel.empty()
        .mix(
            skesa.out.fasta, spades_illumina.out.fasta, 
            spades_iontorrent.out.fasta, medaka.out.fasta
        ).set{ ch_assembly }

    ch_versions = ch_versions.mix(flye.out.versions)
    ch_versions = ch_versions.mix(medaka.out.versions)
    ch_versions = ch_versions.mix(skesa.out.versions)
    ch_versions = ch_versions.mix(spades_illumina.out.versions)
    ch_versions = ch_versions.mix(spades_iontorrent.out.versions)

    emit:
    assembly    = ch_assembly   // channel: [ val(meta), path(fasta) ]
    versions    = ch_versions   // channel: [ versions.yml ]
}
