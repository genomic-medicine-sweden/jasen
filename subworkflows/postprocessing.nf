#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { concatenate_files   } from '../modules/local/jasentool/main.nf'
include { create_yaml         } from '../modules/local/jasentool/main.nf'
include { format_cdm          } from '../modules/local/jasentool/main.nf'
include { export_to_cdm       } from '../modules/local/cdm/main.nf'

workflow CALL_POSTPROCESSING {
    take:
    reference_genome
    reference_genome_idx
    reference_genome_gff
    reference_genome_accession
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
    ch_versions_files

    main:

    ch_versions = Channel.empty()

    concatenate_files(ch_versions_files.collect())

    ch_preprocessing_combined_output
        .join(ch_profiling_combined_output)
        .join(ch_qc_combined_output)
        .join(ch_relatedness_combined_output)
        .join(ch_screening_combined_output)
        .join(ch_typing_combined_output)
        .join(ch_variant_calling_combined_output)
        .set{ ch_combined_output}

    create_yaml(
        ch_combined_output,
        reference_genome,
        reference_genome_idx,
        reference_genome_gff,
        reference_genome_accession,
        tb_grading_rules_bed,
        tbdb_bed,
        concatenate_files.out.concatenated
    )

    format_cdm(create_yaml.out.yaml)

    export_to_cdm(format_cdm.out.json.join(ch_seqrun_meta), species_dir)

    // Fail loudly if any input sample was dropped before producing an analysis YAML.
    ch_preprocessing_combined_output
        .map { it[0] }
        .toList()
        .map { [it] }
        .combine(create_yaml.out.yaml.map { it[0] }.toList().map { [it] })
        .subscribe { expected, actual ->
            def missing = expected - actual
            if (missing) {
                error "Analysis YAML not produced for sample(s): ${missing.sort().join(', ')}"
            }
        }

    emit:
    cdm             = export_to_cdm.out.cdm             // channel: [ path(txt) ]
    yaml            = create_yaml.out.yaml          // channel: [ path(yaml) ]
    versions        = ch_versions                       // channel: [ versions.yml ]
}
