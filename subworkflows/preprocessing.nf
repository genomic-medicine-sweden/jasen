include { assembly_trim_clean       } from '../nextflow-modules/modules/clean/main'
include { save_analysis_metadata    } from '../nextflow-modules/modules/meta/main'

workflow CALL_PREPROCESSING {
    take:
        ch_meta_iontorrent
        ch_meta_illumina

    main:
        // reads trim and clean
        assembly_trim_clean(ch_meta_iontorrent).set{ ch_clean_meta }
        ch_meta_illumina.mix(ch_clean_meta).set{ ch_input_meta }
        ch_input_meta.map { sampleName, reads, platform -> [ sampleName, reads ] }.set{ ch_clean_reads }

        // analysis metadata
        save_analysis_metadata(ch_input_meta)

    emit:
        input_meta  = ch_input_meta                     // channel: [ val(meta), path(meta)]
        reads       = ch_clean_reads                    // channel: [ val(meta), path(json)]
        metadata    = save_analysis_metadata.out.meta   // channel: [ val(meta), path(json)]
}