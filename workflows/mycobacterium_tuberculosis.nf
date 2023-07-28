#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { get_meta                  } from '../methods/get_meta'
include { bracken                   } from '../nextflow-modules/modules/bracken/main'
include { create_analysis_result    } from '../nextflow-modules/modules/prp/main'
include { kraken                    } from '../nextflow-modules/modules/kraken/main'
include { mykrobe                   } from '../nextflow-modules/modules/mykrobe/main'
include { snippy                    } from '../nextflow-modules/modules/snippy/main'
include { tbprofiler                } from '../nextflow-modules/modules/tbprofiler/main'
include { CALL_BACTERIAL_BASE       } from '../workflows/bacterial_base.nf'

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

        mykrobe(CALL_BACTERIAL_BASE.out.reads)

        snippy(CALL_BACTERIAL_BASE.out.reads, genomeReference)

        tbprofiler(CALL_BACTERIAL_BASE.out.reads)
        
        CALL_BACTERIAL_BASE.out.reads.map { sampleName, reads -> [ sampleName, [] ] }.set{ ch_empty }

        CALL_BACTERIAL_BASE.out.quast
            .join(CALL_BACTERIAL_BASE.out.qc)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(CALL_BACTERIAL_BASE.out.metadata)
            .join(mykrobe.out.json)
            .join(snippy.out.vcf)
            .join(tbprofiler.out.json)
            .set{ combinedOutput }

        if ( params.useKraken ) {
            krakenDb = file(params.krakenDb, checkIfExists: true)
            kraken(CALL_BACTERIAL_BASE.out.reads, krakenDb)
            bracken(kraken.out.report, krakenDb).output
            combinedOutput = combinedOutput.join(bracken.out.output)
            //create_analysis_result(combinedOutput)
            ch_versions = ch_versions.mix(kraken.out.versions)
            ch_versions = ch_versions.mix(bracken.out.versions)
        } else {
            emptyBrackenOutput = reads.map { sampleName, reads -> [ sampleName, [] ] }
            combinedOutput = combinedOutput.join(emptyBrackenOutput)
            //create_analysis_result(combinedOutput)
        }

        ch_versions = ch_versions.mix(CALL_BACTERIAL_BASE.out.versions)
        ch_versions = ch_versions.mix(mykrobe.out.versions)
        ch_versions = ch_versions.mix(snippy.out.versions)
        ch_versions = ch_versions.mix(tbprofiler.out.versions)

    emit: 
        //pipeline_result = create_analysis_result.output
        versions        = ch_versions
}