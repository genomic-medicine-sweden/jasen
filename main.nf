#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CALL_STAPHYLOCOCCUS_AUREUS        } from './workflows/staphylococcus_aureus.nf'
include { CALL_MYCOBACTERIUM_TUBERCULOSIS   } from './workflows/mycobacterium_tuberculosis.nf'

workflow {
    if (workflow.profile == "staphylococcus_aureus") {
        CALL_STAPHYLOCOCCUS_AUREUS()
    } else if (workflow.profile == "mycobacterium_tuberculosis") {
        CALL_MYCOBACTERIUM_TUBERCULOSIS()
    }
}
