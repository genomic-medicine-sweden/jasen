#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { get_meta                                  } from '../methods/get_sample_data.nf'
include { get_reads                                 } from '../methods/get_sample_data.nf'
include { annotate_delly                            } from '../nextflow-modules/modules/prp/main.nf'
include { bracken                                   } from '../nextflow-modules/modules/bracken/main.nf'
include { create_analysis_result                    } from '../nextflow-modules/modules/prp/main.nf'
include { add_igv_track as add_tbdb_bed_track       } from '../nextflow-modules/modules/prp/main.nf'
include { add_igv_track as add_grading_bed_track    } from '../nextflow-modules/modules/prp/main.nf'
include { create_cdm_input                          } from '../nextflow-modules/modules/prp/main.nf'
include { create_yaml                               } from '../nextflow-modules/modules/yaml/main.nf'
include { export_to_cdm                             } from '../nextflow-modules/modules/cmd/main.nf'
include { kraken                                    } from '../nextflow-modules/modules/kraken/main.nf'
include { mykrobe                                   } from '../nextflow-modules/modules/mykrobe/main.nf'
include { post_align_qc                             } from '../nextflow-modules/modules/prp/main.nf'
include { snippy                                    } from '../nextflow-modules/modules/snippy/main.nf'
include { tbprofiler as tbprofiler_mergedb          } from '../nextflow-modules/modules/tbprofiler/main.nf'
include { CALL_BACTERIAL_BASE                       } from '../workflows/bacterial_base.nf'

workflow CALL_MYCOBACTERIUM_TUBERCULOSIS {
    // set input data
    inputSamples = file(params.csv, checkIfExists: true)

    // load references 
    referenceGenome = file(params.referenceGenome, checkIfExists: true)
    referenceGenomeDir = file(referenceGenome.getParent(), checkIfExists: true)
    referenceGenomeIdx = file(params.referenceGenomeIdx, checkIfExists: true)
    referenceGenomeGff = file(params.referenceGenomeGff, checkIfExists: true)
    // databases
    coreLociBed = file(params.coreLociBed, checkIfExists: true)
    tbdbBed = file(params.tbdbBed, checkIfExists: true)
    tbdbBedIdx = file(params.tbdbBedIdx, checkIfExists: true)
    tbGradingRulesBed = file(params.tbGradingRulesBed, checkIfExists: true)
    // schemas and values
    targetSampleSize = params.targetSampleSize ? params.targetSampleSize : Channel.value([])

    main:
        ch_versions = Channel.empty()

        CALL_BACTERIAL_BASE( coreLociBed, referenceGenome, referenceGenomeDir, inputSamples, targetSampleSize )

        CALL_BACTERIAL_BASE.out.assembly.set{ch_assembly}
        CALL_BACTERIAL_BASE.out.reads.set{ch_reads}
        CALL_BACTERIAL_BASE.out.quast.set{ch_quast}
        CALL_BACTERIAL_BASE.out.metadata.set{ch_metadata}
        CALL_BACTERIAL_BASE.out.seqplat_meta.set{ch_seqplat_meta}
        CALL_BACTERIAL_BASE.out.seqrun_meta.set{ch_seqrun_meta}
        CALL_BACTERIAL_BASE.out.reads_w_meta.set{ch_input_meta}
        CALL_BACTERIAL_BASE.out.sourmash.set{ch_sourmash}
        CALL_BACTERIAL_BASE.out.ska_build.set{ch_ska}

        mykrobe(ch_reads)

        snippy(ch_reads, referenceGenome)

        tbprofiler_mergedb(ch_reads)

        annotate_delly(tbprofiler_mergedb.out.vcf, tbdbBed, tbdbBedIdx)

        post_align_qc(tbprofiler_mergedb.out.bam, params.referenceGenome, coreLociBed)
        post_align_qc.out.qc.set{ch_qc}

        ch_reads.map { sampleID, reads -> [ sampleID, [] ] }.set{ ch_empty }

        ch_qc_illumina = ch_seqplat_meta.filter { it[1] == 'illumina' } ? ch_qc : ch_empty

        ch_quast
            .join(ch_qc_illumina)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(tbprofiler_mergedb.out.bam)
            .join(tbprofiler_mergedb.out.bai)
            .join(ch_metadata)
            .join(annotate_delly.out.vcf)
            .join(mykrobe.out.csv)
            .join(tbprofiler_mergedb.out.json)
            .set{ combinedOutput }

        if ( params.useKraken ) {
            krakenDb = file(params.krakenDb, checkIfExists: true)
            kraken(ch_reads, krakenDb)
            bracken(kraken.out.report, krakenDb).output
            combinedOutput.join(bracken.out.output).set{ combinedOutput }
            create_analysis_result(combinedOutput, referenceGenome, referenceGenomeIdx, referenceGenomeGff)
            ch_versions = ch_versions.mix(kraken.out.versions)
            ch_versions = ch_versions.mix(bracken.out.versions)
        } else {
            combinedOutput.join(ch_empty).set{ combinedOutput }
            create_analysis_result(combinedOutput, referenceGenome, referenceGenomeIdx, referenceGenomeGff)
        }

        // Add IGV annotation tracks
        add_tbdb_bed_track(create_analysis_result.out.json, params.tbdbBed, params.resistantLociName)
        add_grading_bed_track(add_tbdb_bed_track.out.json, params.tbGradingRulesBed, params.gradingLociName)

        // Create yaml for uploading results to Bonsai
        create_yaml(add_grading_bed_track.out.json.join(ch_sourmash).join(ch_ska), params.speciesDir)

        ch_quast
            .join(ch_qc_illumina)
            .join(ch_empty)
            .set{ cdmInput }

        create_cdm_input(cdmInput)

        export_to_cdm(create_cdm_input.out.json.join(ch_seqrun_meta), params.speciesDir)

        ch_versions = ch_versions.mix(CALL_BACTERIAL_BASE.out.versions)
        ch_versions = ch_versions.mix(create_analysis_result.out.versions)
        ch_versions = ch_versions.mix(mykrobe.out.versions)
        ch_versions = ch_versions.mix(snippy.out.versions)
        ch_versions = ch_versions.mix(tbprofiler_mergedb.out.versions)

    emit: 
        pipeline_result = add_grading_bed_track.out.json
        cdm             = export_to_cdm.out.cdm
        yaml            = create_yaml.out.yaml
        versions        = ch_versions
}
