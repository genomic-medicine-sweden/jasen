#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { create_analysis_result    } from '../modules/local/prp/main.nf'
include { create_cdm_input          } from '../modules/local/prp/main.nf'
include { create_prp_yaml           } from '../modules/local/yaml/prp/main.nf'
include { export_to_cdm             } from '../modules/local/cdm/main.nf'

workflow CALL_POSTPROCESSING {
    take:
    reference_genome
    reference_genome_idx
    reference_genome_gff
    species_dir
    ch_chewbbaca
    ch_post_align_qc
    ch_preprocessing_combined_output
    ch_profiling_combined_output
    ch_qc_combined_output
    ch_quast
    ch_relatedness_combined_output
    ch_screening_combined_output
    ch_seqrun_meta
    ch_typing_combined_output
    ch_variant_calling_combined_output

    main:

    ch_versions = Channel.empty()

    // POSTPROCESSING
    ch_preprocessing_combined_output
        .join(ch_profiling_combined_output)
        .join(ch_qc_combined_output)
        .join(ch_relatedness_combined_output)
        .join(ch_screening_combined_output)
        .join(ch_typing_combined_output)
        .join(ch_variant_calling_combined_output)
        .set{ ch_combined_output}

    create_prp_yaml(
        ch_combined_output,
        reference_genome,
        reference_genome_idx,
        reference_genome_gff
    )

    create_analysis_result(create_prp_yaml.out.yaml)

    ch_quast
        .join(ch_post_align_qc)
        .join(ch_chewbbaca)
        .set{ ch_cdm_input }

    create_cdm_input(ch_cdm_input)

    export_to_cdm(create_cdm_input.out.json.join(ch_seqrun_meta), species_dir)

    emit:
    pipeline_result = create_analysis_result.out.json   // channel: [ path(json) ]
    cdm             = export_to_cdm.out.cdm             // channel: [ path(txt) ]
    yaml            = create_prp_yaml.out.yaml          // channel: [ path(yaml) ]
    versions        = ch_versions                       // channel: [ versions.yml ]
}
