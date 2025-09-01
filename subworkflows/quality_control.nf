#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { bracken                                    } from '../modules/nf-core/bracken/main.nf'
include { bwa_mem as bwa_mem_ref                     } from '../modules/nf-core/bwa/main.nf'
include { fastqc                                     } from '../modules/nf-core/fastqc/main.nf'
include { gambitcore                                 } from '../modules/local/gambitcore/main.nf'
include { kraken                                     } from '../modules/nf-core/kraken/main.nf'
include { minimap2_align as minimap2_align_ref       } from '../modules/nf-core/minimap2/main.nf'       
include { nanoplot                                   } from '../modules/nf-core/nanoplot/main.nf'
include { post_align_qc                              } from '../modules/local/prp/main.nf'
include { quast                                      } from '../modules/nf-core/quast/main.nf'
include { samtools_coverage as samtools_coverage_ref } from '../modules/nf-core/samtools/main.nf'
include { samtools_index as samtools_index_ref       } from '../modules/nf-core/samtools/main.nf'
include { samtools_sort as samtools_sort_ref         } from '../modules/nf-core/samtools/main.nf'

workflow CALL_QUALITY_CONTROL {
    take:
    core_loci_bed
    gambit_db
    kraken_db
    reference_genome
    reference_genome_dir
    reference_genome_idx
    ch_assembly
    ch_empty
    ch_reads

    main:

    ch_versions = Channel.empty()

    // evaluate assembly completeness
    gambitcore(ch_assembly, gambit_db)

    // evaluate assembly quality 
    quast(ch_assembly, reference_genome)

    // qc processing - short read
    fastqc(ch_reads)

    bwa_mem_ref(ch_reads, reference_genome_dir)

    // qc processing - long read
    nanoplot(ch_reads)

    minimap2_align_ref(ch_reads, reference_genome_idx)
    samtools_sort_ref(minimap2_align_ref.out.sam)

    if (params.reference_genome) {
        samtools_sort_ref.out.bam.mix(bwa_mem_ref.out.bam).set{ ch_ref_bam }
        samtools_index_ref(ch_ref_bam).bai.set{ ch_ref_bai }
        post_align_qc(ch_ref_bam, reference_genome, core_loci_bed).json.set{ ch_post_align_qc }
        samtools_coverage_ref(samtools_sort_ref.out.bam).txt.set{ ch_samtools_cov_ref }
        ch_versions = ch_versions.mix(samtools_coverage_ref.out.versions)
        ch_versions = ch_versions.mix(samtools_index_ref.out.versions)
        ch_versions = ch_versions.mix(samtools_sort_ref.out.versions)
    } else {
        ch_empty.set{ ch_ref_bam }
        ch_empty.set{ ch_ref_bai }
        ch_empty.set{ ch_post_align_qc }
        ch_empty.set{ ch_samtools_cov_ref }
    }

    if ( params.use_kraken ) {
        kraken(ch_reads, kraken_db)
        bracken(kraken.out.report, kraken_db)
        bracken.out.output.set{ ch_kraken }
        ch_versions = ch_versions.mix(kraken.out.versions)
        ch_versions = ch_versions.mix(bracken.out.versions)
    } else {
        ch_empty.set{ ch_kraken }
    }

    ch_ref_bam
        .join(ch_ref_bai)
        .join(gambitcore.out.tsv)
        .join(ch_kraken)
        .join(ch_post_align_qc)
        .join(quast.out.tsv)
        .set { ch_combined_output }

    ch_versions = ch_versions.mix(bwa_mem_ref.out.versions)
    ch_versions = ch_versions.mix(fastqc.out.versions)
    ch_versions = ch_versions.mix(minimap2_align_ref.out.versions)
    ch_versions = ch_versions.mix(nanoplot.out.versions)

    emit:
    combined_output     = ch_combined_output            // channel: [ val(meta), val(bam), val(bai), path(txt), path(json), path(tsv) ]
    bam                 = ch_ref_bam                    // channel: [ val(meta), path(bam) ]
    bai                 = ch_ref_bai                    // channel: [ val(meta), path(bai) ]
    fastqc              = fastqc.out.output             // channel: [ val(meta), path(txt) ]
    gambitcore          = gambitcore.out.tsv            // channel: [ val(meta), path(tsv) ]
    kraken              = ch_kraken                     // channel: [ val(meta), path(fasta) ]
    nanoplot_html       = nanoplot.out.html             // channel: [ val(meta), path(html) ]
    nanoplot_txt        = nanoplot.out.txt              // channel: [ val(meta), path(txt) ]
    post_align_qc       = ch_post_align_qc              // channel: [ val(meta), path(fasta) ]
    quast               = quast.out.tsv                 // channel: [ val(meta), path(tsv) ]
    samtools_cov_ref    = ch_samtools_cov_ref           // channel: [ val(meta), path(txt) ]
    versions            = ch_versions                   // channel: [ versions.yml ]
}
