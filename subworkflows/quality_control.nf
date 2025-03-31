#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { bracken                                    } from '../modules/nf-core/bracken/main.nf'
include { bwa_mem as bwa_mem_ref                     } from '../modules/nf-core/bwa/main.nf'
include { fastqc                                     } from '../modules/nf-core/fastqc/main.nf'
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
    kraken_db
    reference_genome
    reference_genome_dir
    reference_genome_idx
    ch_assembly
    ch_reads
    ch_reads_w_meta
    ch_versions

    main:
    // evaluate assembly quality 
    quast(ch_assembly, reference_genome)

    // qc processing - short read
    fastqc(ch_reads_w_meta)
    bwa_mem_ref(ch_reads_w_meta, reference_genome_dir)

    post_align_qc(bwa_mem_ref.out.bam, reference_genome, core_loci_bed)

    // qc processing - long read
    nanoplot(ch_reads_w_meta)

    minimap2_align_ref(ch_reads_w_meta, reference_genome_idx)
    samtools_sort_ref(minimap2_align_ref.out.sam)
    samtools_coverage_ref(samtools_sort_ref.out.bam)

    samtools_sort_ref.out.bam.mix(bwa_mem_ref.out.bam).set{ ch_ref_bam }

    samtools_index_ref(ch_ref_bam)

    if ( params.use_kraken ) {
        kraken(ch_reads, kraken_db)
        bracken(kraken.out.report, kraken_db)
        bracken.out.output.set{ ch_kraken }
        ch_versions = ch_versions.mix(kraken.out.versions)
        ch_versions = ch_versions.mix(bracken.out.versions)
    } else {
        ch_empty.set{ ch_kraken }
    }

    ch_versions = ch_versions.mix(bwa_mem_ref.out.versions)
    ch_versions = ch_versions.mix(fastqc.out.versions)
    ch_versions = ch_versions.mix(minimap2_align_ref.out.versions)
    ch_versions = ch_versions.mix(nanoplot.out.versions)
    ch_versions = ch_versions.mix(post_align_qc.out.versions)
    ch_versions = ch_versions.mix(samtools_coverage_ref.out.versions)
    ch_versions = ch_versions.mix(samtools_index_ref.out.versions)
    ch_versions = ch_versions.mix(samtools_sort_ref.out.versions)

    emit:
    bam                 = ch_ref_bam                    // channel: [ val(meta), path(bam) ]
    bai                 = samtools_index_ref.out.bai    // channel: [ val(meta), path(bai) ]
    quast               = quast.out.tsv                 // channel: [ val(meta), path(tsv)]
    kraken              = ch_kraken                     // channel: [ val(meta), path(fasta)]
    post_align_qc       = post_align_qc.out.json        // channel: [ val(meta), path(fasta)]
    nanoplot_html       = nanoplot.out.html             // channel: [ val(meta), path(html) ]
    samtools_cov_txt    = samtools_coverage_ref.out.txt // channel: [ val(meta), path(txt)]
    versions            = ch_versions                   // channel: [ versions.yml ]
}
