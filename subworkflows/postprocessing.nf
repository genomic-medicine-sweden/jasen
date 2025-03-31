#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { create_analysis_result    } from '../modules/local/prp/main.nf'
include { export_to_cdm             } from '../modules/local/cmd/main.nf'
include { create_prp_yaml           } from '../modules/local/yaml/prp/main.nf'

workflow CALL_POSTPROCESSING {
    take:
    reference_genome
    reference_genome_idx
    reference_genome_gff
    ch_post_align_qc
    ch_preprocessing_combined_output
    ch_qc_combined_output
    ch_quast
    ch_relatedness_combined_output
    ch_screening_combined_output
    ch_seqrun_meta
    ch_typing_combined_output
    ch_variant_calling_combined_output
    ch_versions

    main:
    // POSTPROCESSING
    ch_preprocessing_combined_output
        .join(ch_qc_combined_output)
        .join(ch_typing_combined_output)
        .join(ch_screening_combined_output)
        .join(ch_relatedness_combined_output)
        .join(ch_variant_calling_combined_output)
        .set{ combined_output }

    create_prp_yaml(combined_output, reference_genome, reference_genome_idx, reference_genome_gff)

    create_analysis_result(create_prp_yaml.out.yaml)

    ch_quast
        .join(ch_post_align_qc)
        .join(ch_chewbbaca)
        .set{ ch_cdm_input }

    create_cdm_input(ch_cdm_input)

    export_to_cdm(create_cdm_input.out.json.join(ch_seqrun_meta), params.species_dir)

    emit:
    pipeline_result = create_analysis_result.output // channel: [ path(json) ]
    cdm             = export_to_cdm.output          // channel: [ path(txt) ]
    versions        = ch_versions                   // channel: [ versions.yml ]
}
