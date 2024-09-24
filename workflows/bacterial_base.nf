#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { get_seqrun_meta                       } from '../methods/get_seqrun_meta.nf'
include { assembly_trim_clean                   } from '../nextflow-modules/modules/clean/main.nf'
include { bwa_mem as bwa_mem_ref                } from '../nextflow-modules/modules/bwa/main.nf'
include { flye                                  } from '../nextflow-modules/modules/flye/main.nf'
include { medaka                                } from '../nextflow-modules/modules/medaka/main.nf'
include { post_align_qc                         } from '../nextflow-modules/modules/prp/main.nf'
include { quast                                 } from '../nextflow-modules/modules/quast/main.nf'
include { samtools_index as samtools_index_ref  } from '../nextflow-modules/modules/samtools/main.nf'
include { save_analysis_metadata                } from '../nextflow-modules/modules/meta/main.nf'
include { skesa                                 } from '../nextflow-modules/modules/skesa/main.nf'
include { sourmash                              } from '../nextflow-modules/modules/sourmash/main.nf'
include { spades_illumina                       } from '../nextflow-modules/modules/spades/main.nf'
include { spades_iontorrent                     } from '../nextflow-modules/modules/spades/main.nf'

workflow CALL_BACTERIAL_BASE {
    take:
        coreLociBed
        referenceGenome
        referenceGenomeDir
        ch_meta_iontorrent
        ch_meta_illumina
        ch_meta_nanopore
    
    main:
        ch_versions = Channel.empty()

        // reads trim and clean
        assembly_trim_clean(ch_meta_iontorrent).set{ ch_clean_meta }
        ch_meta_illumina.mix(ch_clean_meta, ch_meta_nanopore).set{ ch_input_meta }
        ch_input_meta.map { sampleID, reads, platform -> [ sampleID, reads ] }.set{ ch_reads }

        if ( params.cronCopy || params.devMode) {
            Channel.fromPath(params.csv).splitCsv(header:true)
                .map{ row -> get_seqrun_meta(row) }
                .set{ ch_seqrun_meta }
        } else {
            ch_reads.map { sampleID, reads -> [ sampleID, [], [], [] ] }.set{ ch_seqrun_meta }
        }

        // analysis metadata
        save_analysis_metadata(ch_input_meta.join(ch_seqrun_meta))

        // assembly
        skesa(ch_input_meta)
        spades_illumina(ch_input_meta)
        spades_iontorrent(ch_input_meta)
        flye(ch_input_meta)
        medaka(ch_input_meta, flye.out.fasta)

        Channel.empty().mix(skesa.out.fasta, spades_illumina.out.fasta, spades_iontorrent.out.fasta, medaka.out.fasta).set{ ch_assembly }

        // evaluate assembly quality 
        quast(ch_assembly, referenceGenome)

        // qc processing
        bwa_mem_ref(ch_reads, referenceGenomeDir)
        samtools_index_ref(bwa_mem_ref.out.bam)

        post_align_qc(bwa_mem_ref.out.bam, referenceGenome, coreLociBed)

        sourmash(ch_assembly)

        ch_versions = ch_versions.mix(bwa_mem_ref.out.versions)
        ch_versions = ch_versions.mix(flye.out.versions)
        ch_versions = ch_versions.mix(medaka.out.versions)
        ch_versions = ch_versions.mix(quast.out.versions)
        ch_versions = ch_versions.mix(samtools_index_ref.out.versions)
        ch_versions = ch_versions.mix(skesa.out.versions)
        ch_versions = ch_versions.mix(sourmash.out.versions)
        ch_versions = ch_versions.mix(spades_illumina.out.versions)
        ch_versions = ch_versions.mix(spades_iontorrent.out.versions)

    emit:
        assembly    = ch_assembly                       // channel: [ val(meta), path(fasta)]
        input_meta  = ch_input_meta                     // channel: [ val(meta), path(meta)]
        bam         = bwa_mem_ref.out.bam               // channel: [ val(meta), path(bam)]
        bai         = samtools_index_ref.out.bai        // channel: [ val(meta), path(bai)]
        metadata    = save_analysis_metadata.out.meta   // channel: [ val(meta), path(json)]
        qc          = post_align_qc.out.qc              // channel: [ val(meta), path(fasta)]
        quast       = quast.out.qc                      // channel: [ val(meta), path(qc)]
        reads       = ch_reads                          // channel: [ val(meta), path(json)]
        seqrun_meta = ch_seqrun_meta                    // channel: [ val(meta), val(json), val(json)]
        sourmash    = sourmash.out.signature            // channel: [ val(meta), path(signature)]
        versions    = ch_versions                       // channel: [ versions.yml ]
}
