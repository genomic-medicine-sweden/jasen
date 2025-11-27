#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { get_reads                 } from '../methods/get_sample_data.nf'
include { get_seqrun_meta           } from '../methods/get_seqrun_meta.nf'
include { assembly_trim_clean       } from '../modules/local/clean/main.nf'
include { hostile                   } from '../modules/nf-core/hostile/main.nf'
include { save_analysis_metadata    } from '../modules/local/meta/main.nf'
include { seqtk_sample              } from '../modules/nf-core/seqtk/main.nf'

workflow CALL_PREPROCESSING {
    take:
    assay
    hostile_dir
    hostile_idx
    input_samples
    platform
    release_life_cycle

    main:

    ch_versions = Channel.empty()

    // PREPROCESSING
    // Create channel for reads
    Channel.fromPath(input_samples)
        .splitCsv(header:true)
        .tap{ ch_raw_input }
        .map{ row -> get_reads(row) }
        .set{ ch_raw_reads }

    if ( params.use_hostile ) {
        // remove human reads
        hostile( ch_raw_reads, hostile_dir, hostile_idx ).reads.set{ ch_depleted_reads }
        ch_versions = ch_versions.mix(hostile.out.versions)
    } else {
        ch_raw_reads.set{ ch_depleted_reads }
    }

    if ( params.target_sample_size ) {
        // downsample reads
        seqtk_sample( ch_depleted_reads, target_sample_size ).reads.set{ ch_depleted_sampled_reads }
        ch_versions = ch_versions.mix(seqtk_sample.out.versions)
    } else {
        ch_depleted_reads.set{ ch_depleted_sampled_reads }
    }

    // reads trim and clean and recreate reads channel if the reads were filtered or downsampled
    assembly_trim_clean(ch_depleted_sampled_reads).set { ch_clean_reads }
    Channel.empty()
        .mix( ch_depleted_sampled_reads, ch_clean_reads )  // if samples are filtered or downsampled
        .set{ ch_reads }                                   // create reads channel

    Channel.fromPath(input_samples).splitCsv(header:true)
        .map{ row -> get_seqrun_meta(row) }
        .tap{ ch_seqrun_meta }
        .map{ id, sequencing_run, lims_id, sample_name -> [ id, lims_id, sample_name ]}
        .set{ ch_id_meta }

    // create empty channel containing only sample_id
    ch_reads.map{ sample_id, reads -> [ sample_id, [] ] }.set{ ch_sample_id }

    // analysis metadata
    save_analysis_metadata(ch_reads.join(ch_seqrun_meta), assay, platform, release_life_cycle)

    ch_id_meta.join(save_analysis_metadata.out.json).set{ ch_combined_output }

    emit:
    combined_output     = ch_combined_output                // channel: [ val(meta), val(meta), val(meta), path(json) ]
    sample_id           = ch_sample_id                      // channel: [ val(meta) ]
    id_meta             = ch_id_meta                        // channel: [ val(meta), val(meta), val(meta), val(meta) ]
    nextflow_run_info   = save_analysis_metadata.out.json   // channel: [ val(meta), path(json) ]
    reads               = ch_reads                          // channel: [ val(meta), path(json) ]
    seqrun_meta         = ch_seqrun_meta                    // channel: [ val(meta), val(json), val(json) ]
    versions            = ch_versions                       // channel: [ versions.yml ]
}
