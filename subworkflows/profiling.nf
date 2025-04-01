#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { annotate_delly                            } from '../modules/local/prp/main.nf'
include { mykrobe                                   } from '../modules/nf-core/mykrobe/main.nf'
include { snippy                                    } from '../modules/nf-core/snippy/main.nf'
include { tbprofiler as tbprofiler_mergedb          } from '../modules/nf-core/tbprofiler/main.nf'

workflow CALL_PROFILING {
    take:
    core_loci_bed
    reference_genome
    tbdb_bed
    tbdb_bed_idx
    ch_assembly
    ch_reads

    main:

    ch_versions = Channel.empty()

    // PROFILING
    mykrobe(ch_reads)

    snippy(ch_reads, reference_genome)

    tbprofiler_mergedb(ch_reads)

    annotate_delly(tbprofiler_mergedb.out.vcf, tbdb_bed, tbdb_bed_idx)

    mykrobe.out.csv
        .join(tbprofiler_mergedb.out.json)
        .set{ ch_combined_output }

    ch_versions = ch_versions.mix(annotate_delly.out.versions)
    ch_versions = ch_versions.mix(mykrobe.out.versions)
    ch_versions = ch_versions.mix(snippy.out.versions)
    ch_versions = ch_versions.mix(tbprofiler_mergedb.out.versions)

    emit:
    bam             = tbprofiler_mergedb.out.bam    // channel: [ val(meta), path(bam) ]
    bai             = tbprofiler_mergedb.out.bai    // channel: [ val(meta), path(bai) ]
    combined_output = ch_combined_output            // channel: [ val(meta), path(csv), path(json) ]
    mykrobe         = mykrobe.out.csv               // channel: [ val(meta), path(csv) ]
    tbprofiler      = tbprofiler_mergedb.out.json   // channel: [ val(meta), path(json) ]
    vcf             = annotate_delly.out.vcf        // channel: [ val(meta), path(json) ]
    versions        = ch_versions                   // channel: [ versions.yml ]
}
