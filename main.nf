#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CALL_STAPHYLOCOCCUS_AUREUS        } from './workflows/staphylococcus_aureus.nf'
include { CALL_MYCOBACTERIUM_TUBERCULOSIS   } from './workflows/mycobacterium_tuberculosis.nf'
include { CALL_KLEBSIELLA_PNEUMONIAE        } from './workflows/klebsiella_pneumoniae.nf'
include { CALL_ESCHERICHIA_COLI             } from './workflows/escherichia_coli.nf'
include { CALL_STREPTOCOCCUS_PYOGENES       } from './workflows/streptococcus_pyogenes.nf'
include { CALL_STREPTOCOCCUS                } from './workflows/streptococcus.nf'

workflow {
    if (params.species == "staphylococcus aureus") {
        CALL_STAPHYLOCOCCUS_AUREUS()
    } else if (params.species == "mycobacterium tuberculosis") {
        CALL_MYCOBACTERIUM_TUBERCULOSIS()
    } else if (params.species == "klebesiella pneumoniae") {
        CALL_KLEBSIELLA_PNEUMONIAE()
    } else if (params.species == "escherichia coli") {
        CALL_ESCHERICHIA_COLI()
    } else if (params.species == "streptococcus pyogenes") {
        CALL_STREPTOCOCCUS_PYOGENES()
    } else if (params.species == "streptococcus") {
        CALL_STREPTOCOCCUS()
    }
}
