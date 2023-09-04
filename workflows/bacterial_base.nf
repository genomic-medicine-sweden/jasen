#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { get_meta                              } from '../methods/get_meta'
include { quast                                 } from '../nextflow-modules/modules/quast/main'
include { skesa                                 } from '../nextflow-modules/modules/skesa/main'
include { spades_illumina                       } from '../nextflow-modules/modules/spades/main'
include { spades_iontorrent                     } from '../nextflow-modules/modules/spades/main'
include { bwa_mem as bwa_mem_ref                } from '../nextflow-modules/modules/bwa/main'
include { samtools_index as samtools_index_ref  } from '../nextflow-modules/modules/samtools/main'
include { post_align_qc                         } from '../nextflow-modules/modules/qc/main'
include { assembly_trim_clean                   } from '../nextflow-modules/modules/clean/main'
include { save_analysis_metadata                } from '../nextflow-modules/modules/meta/main'
include { sourmash                              } from '../nextflow-modules/modules/sourmash/main'

workflow CALL_BACTERIAL_BASE {
    take:
        coreLociBed
        genomeReference
        genomeReferenceDir
        ch_meta_iontorrent
        ch_meta_illumina
    
    main:
        ch_versions = Channel.empty()

        // reads trim and clean
        assembly_trim_clean(ch_meta_iontorrent).set{ ch_clean_meta }
        ch_meta_illumina.mix(ch_clean_meta).set{ ch_input_meta }
        ch_input_meta.map { sampleName, reads, platform -> [ sampleName, reads ] }.set{ ch_reads }

        // analysis metadata
        save_analysis_metadata(ch_input_meta)

        // assembly
        skesa(ch_input_meta)
        spades_illumina(ch_input_meta)
        spades_iontorrent(ch_input_meta)

        Channel.empty().mix(skesa.out.fasta, spades_illumina.out.fasta, spades_iontorrent.out.fasta).set{ ch_assembly }

        // evaluate assembly quality 
        quast(ch_assembly, genomeReference)

        // qc processing
        bwa_mem_ref(ch_reads, genomeReferenceDir)
        samtools_index_ref(bwa_mem_ref.out.bam)

        bwa_mem_ref.out.bam
            .join(samtools_index_ref.out.bai)
            .multiMap { id, bam, bai -> 
                bam: tuple(id, bam)
                bai: bai
            }
            .set{ post_align_qc_ch }

        post_align_qc(post_align_qc_ch.bam, post_align_qc_ch.bai, coreLociBed)

        sourmash(ch_assembly)

        ch_versions = ch_versions.mix(bwa_mem_ref.out.versions)
        ch_versions = ch_versions.mix(quast.out.versions)
        ch_versions = ch_versions.mix(samtools_index_ref.out.versions)
        ch_versions = ch_versions.mix(skesa.out.versions)
        ch_versions = ch_versions.mix(sourmash.out.versions)
        ch_versions = ch_versions.mix(spades_illumina.out.versions)
        ch_versions = ch_versions.mix(spades_iontorrent.out.versions)

    emit:
        assembly    = ch_assembly                       // channel: [ val(meta), path(fasta)]
        input_meta  = ch_input_meta                     // channel: [ val(meta), path(meta)]
        metadata    = save_analysis_metadata.out.meta   // channel: [ val(meta), path(json)]
        qc          = post_align_qc.out.qc              // channel: [ val(meta), path(fasta)]
        quast       = quast.out.qc                      // channel: [ val(meta), path(qc)]
        reads       = ch_reads                          // channel: [ val(meta), path(json)]
        sourmash    = sourmash.out.signature            // channel: [ val(meta), path(signature)]
        versions    = ch_versions                       // channel: [ versions.yml ]
}