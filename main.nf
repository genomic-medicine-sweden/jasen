#!/usr/bin/env nextflow

nextflow.enable.dsl=2

if (!nextflow.version.matches('>=25.04.0')) {
    nextflow.preview.topic = true
}

include { CALL_MYCOBACTERIUM_TUBERCULOSIS   } from './workflows/mycobacterium_tuberculosis.nf'
include { CALL_BACTERIAL_GENERAL            } from './workflows/bacterial_general.nf'

workflow {
    if (params.species == "mycobacterium tuberculosis") {
        CALL_MYCOBACTERIUM_TUBERCULOSIS()
    } else {
        CALL_BACTERIAL_GENERAL()
    } 
}
