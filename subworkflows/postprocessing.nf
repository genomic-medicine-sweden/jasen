include { bracken                   } from '../nextflow-modules/modules/bracken/main'
include { create_analysis_result    } from '../nextflow-modules/modules/prp/main'
include { export_to_cdm             } from '../nextflow-modules/modules/cmd/main'
include { kraken                    } from '../nextflow-modules/modules/kraken/main'

workflow CALL_POSTPROCESSING {
    take:
        ch_chewbbaca_split_results  // channel: [ val(meta), val(ch_post_align_qc) ]
        ch_combinedOutput          // channel: [ val(meta), val(ch_post_align_qc) ]
        ch_post_align_qc            // channel: [ val(meta), val(ch_post_align_qc) ]
        ch_quast                    // channel: [ val(meta), val(quast) ]
        ch_reads                    // channel: [ val(meta), val(reads) ]

    main:
        ch_versions = Channel.empty()

        // end point
        export_to_cdm(ch_chewbbaca_split_results.join(ch_quast).join(ch_post_align_qc))

        // Using kraken for species identificaiton
        if ( params.useKraken ) {
            krakenDb = file(params.krakenDb, checkIfExists: true)
            kraken(ch_reads, krakenDb)
            bracken(kraken.out.report, krakenDb).output
            ch_combinedOutput = ch_combinedOutput.join(bracken.out.output)
            create_analysis_result(ch_combinedOutput)
            ch_versions = ch_versions.mix(kraken.out.versions)
            ch_versions = ch_versions.mix(bracken.out.versions)
        } else {
            emptyBrackenOutput = reads.map { sampleName, reads -> [ sampleName, [] ] }
            ch_combinedOutput = ch_combinedOutput.join(emptyBrackenOutput)
            create_analysis_result(ch_combinedOutput)
        }

    emit:
        pipeline_result = create_analysis_result.output // channel: [ path(json) ]
        cdm             = export_to_cdm.output          // channel: [ path(txt) ]
        versions        = ch_versions                   // channel: [ versions.yml ]
}