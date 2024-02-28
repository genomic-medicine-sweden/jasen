#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { get_meta                          } from '../methods/get_meta.nf'
include { get_seqrun_meta                   } from '../methods/get_seqrun_meta.nf'
include { bracken                           } from '../nextflow-modules/modules/bracken/main.nf'
include { copy_to_cron                      } from '../nextflow-modules/modules/cron/main.nf'
include { create_analysis_result            } from '../nextflow-modules/modules/prp/main.nf'
include { create_cdm_input                  } from '../nextflow-modules/modules/prp/main.nf'
include { export_to_cdm                     } from '../nextflow-modules/modules/cmd/main.nf'
include { kraken                            } from '../nextflow-modules/modules/kraken/main.nf'
include { mykrobe                           } from '../nextflow-modules/modules/mykrobe/main.nf'
include { snippy                            } from '../nextflow-modules/modules/snippy/main.nf'
include { tbprofiler as tbprofiler_tbdb     } from '../nextflow-modules/modules/tbprofiler/main.nf'
include { tbprofiler as tbprofiler_mergedb  } from '../nextflow-modules/modules/tbprofiler/main.nf'
include { CALL_BACTERIAL_BASE               } from '../workflows/bacterial_base.nf'

workflow CALL_MYCOBACTERIUM_TUBERCULOSIS {
    Channel.fromPath(params.csv).splitCsv(header:true)
        .map{ row -> get_meta(row) }
        .branch {
        iontorrent: it[2] == "iontorrent"
        illumina: it[2] == "illumina"
        }
        .set{ ch_meta }

    // load references 
    genomeReference = file(params.genomeReference, checkIfExists: true)
    genomeReferenceDir = file(genomeReference.getParent(), checkIfExists: true)
    // databases
    coreLociBed = file(params.coreLociBed, checkIfExists: true)

    main:
        ch_versions = Channel.empty()

        CALL_BACTERIAL_BASE( coreLociBed, genomeReference, genomeReferenceDir, ch_meta.iontorrent, ch_meta.illumina )

        CALL_BACTERIAL_BASE.out.assembly.set{ch_assembly}
        CALL_BACTERIAL_BASE.out.reads.set{ch_reads}
        CALL_BACTERIAL_BASE.out.quast.set{ch_quast}
        CALL_BACTERIAL_BASE.out.qc.set{ch_qc}
        CALL_BACTERIAL_BASE.out.metadata.set{ch_metadata}
        CALL_BACTERIAL_BASE.out.input_meta.set{ch_input_meta}

        mykrobe(ch_reads)

        snippy(ch_reads, genomeReference)

        tbprofiler_tbdb(ch_reads)
        tbprofiler_mergedb(ch_reads)

        ch_reads.map { sampleName, reads -> [ sampleName, [] ] }.set{ ch_empty }

        ch_quast
            .join(ch_qc)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_metadata)
            .join(mykrobe.out.csv)
            .join(tbprofiler_tbdb.out.json)
            .set{ combinedOutput }

        if ( params.useKraken ) {
            krakenDb = file(params.krakenDb, checkIfExists: true)
            kraken(ch_reads, krakenDb)
            bracken(kraken.out.report, krakenDb).output
            combinedOutput.join(bracken.out.output).set{ combinedOutput }
            create_analysis_result(combinedOutput)
            ch_versions = ch_versions.mix(kraken.out.versions)
            ch_versions = ch_versions.mix(bracken.out.versions)
        } else {
            combinedOutput.join(ch_empty).set{ combinedOutput }
            create_analysis_result(combinedOutput)
        }

        ch_quast
            .join(ch_qc)
            .join(ch_empty)
            .set{ cdmOutput }

        create_cdm_input(cdmOutput)

        if ( params.cronCopy ) {
            Channel.fromPath(params.csv).splitCsv(header:true)
                .map{ row -> get_seqrun_meta(row) }
                .set{ ch_seqrun_meta }
        } else {
            ch_empty.join(ch_empty).set{ ch_seqrun_meta }
        }

        export_to_cdm(create_cdm_input.out.json.join(ch_seqrun_meta), params.speciesDir)

        copy_to_cron(create_analysis_result.out.json.join(export_to_cdm.out.cdm))

        ch_versions = ch_versions.mix(CALL_BACTERIAL_BASE.out.versions)
        ch_versions = ch_versions.mix(mykrobe.out.versions)
        ch_versions = ch_versions.mix(snippy.out.versions)
        ch_versions = ch_versions.mix(tbprofiler_tbdb.out.versions)

    emit: 
        pipeline_result = create_analysis_result.out.json
        cdm             = export_to_cdm.out.cdm
        cron_json       = copy_to_cron.out.json
        cron_cdm        = copy_to_cron.out.cdm
        versions        = ch_versions
}