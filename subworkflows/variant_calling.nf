#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { bwa_index                                 } from '../modules/nf-core/bwa/main.nf'
include { bwa_mem as bwa_mem_assembly               } from '../modules/nf-core/bwa/main.nf'
include { freebayes                                 } from '../modules/nf-core/freebayes/main.nf'   
include { samtools_index as samtools_index_assembly } from '../modules/nf-core/samtools/main.nf'

workflow CALL_VARIANT_CALLING {
    take:
    ch_assembly
    ch_reads

    main:

    ch_versions = Channel.empty()

    // VARIANT CALLING
    bwa_index(ch_assembly)

    // create input map channels for bwa on assembly
    ch_reads
        .join(bwa_index.out.index)
        .multiMap { id, reads, index -> 
            reads: tuple(id, reads)
            index: index
        }
        .set{ ch_bwa_mem_assembly_map }
    bwa_mem_assembly(ch_bwa_mem_assembly_map.reads, ch_bwa_mem_assembly_map.index)

    samtools_index_assembly(bwa_mem_assembly.out.bam)

    // construct freebayes input channels
    ch_bam_bai = bwa_mem_assembly.out.bam
        .join(samtools_index_assembly.out.bai)
        .set{ ch_bam_bai }

    freebayes(ch_assembly.join(ch_bam_bai))

    ch_versions = ch_versions.mix(bwa_index.out.versions)
    ch_versions = ch_versions.mix(bwa_mem_assembly.out.versions)
    ch_versions = ch_versions.mix(freebayes.out.versions)
    ch_versions = ch_versions.mix(samtools_index_assembly.out.versions)

    emit:
    vcf         = freebayes.out.vcf     // channel: [ val(meta), path(vcf) ]
    versions    = ch_versions           // channel: [ versions.yml ]
}
