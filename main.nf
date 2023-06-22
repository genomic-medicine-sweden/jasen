#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CALL_BACTERIAL_DEFAULT  } from './workflows/bacterial_default.nf'

workflow {
    CALL_BACTERIAL_DEFAULT()  
}
