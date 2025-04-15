#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { format_jasen      } from '../modules/local/prp/main.nf'
include { format_cdm        } from '../modules/local/prp/main.nf'
include { create_prp_yaml   } from '../modules/local/yaml/prp/main.nf'
include { export_to_cdm     } from '../modules/local/cdm/main.nf'

workflow CALL_POSTPROCESSING {
    take:
    reference_genome
    reference_genome_idx
    reference_genome_gff
    species_dir
    tb_grading_rules_bed
    tbdb_bed
    ch_preprocessing_combined_output
    ch_profiling_combined_output
    ch_qc_combined_output
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
        reference_genome_gff,
        tb_grading_rules_bed,
        tbdb_bed
    )

    format_jasen(create_prp_yaml.out.yaml)

    format_cdm(create_prp_yaml.out.yaml)

    export_to_cdm(format_cdm.out.json.join(ch_seqrun_meta), species_dir)

    emit:
    pipeline_result = format_jasen.out.json             // channel: [ path(json) ]
    cdm             = export_to_cdm.out.cdm             // channel: [ path(txt) ]
    yaml            = create_prp_yaml.out.yaml          // channel: [ path(yaml) ]
    versions        = ch_versions                       // channel: [ versions.yml ]
}
