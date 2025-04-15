#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CALL_MYCOBACTERIUM_TUBERCULOSIS   } from './workflows/mycobacterium_tuberculosis.nf'
include { CALL_BACTERIAL_GENERAL            } from './workflows/bacterial_general.nf'

workflow {
    if (params.species == "mycobacterium tuberculosis") {
        CALL_MYCOBACTERIUM_TUBERCULOSIS()
    } else {
        CALL_BACTERIAL_GENERAL()
    } 
}
