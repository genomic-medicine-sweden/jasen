#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { bwa_index                                 } from '../modules/nf-core/bwa/main.nf'
include { bwa_mem as bwa_mem_assembly               } from '../modules/nf-core/bwa/main.nf'
include { freebayes                                 } from '../modules/nf-core/freebayes/main.nf'
include { minimap2_align as minimap2_align_assembly } from '../modules/nf-core/minimap2/main.nf'       
include { minimap2_index                            } from '../modules/nf-core/minimap2/main.nf'       
include { samtools_index as samtools_index_assembly } from '../modules/nf-core/samtools/main.nf'
include { samtools_sort as samtools_sort_assembly   } from '../modules/nf-core/samtools/main.nf'

workflow CALL_VARIANT_CALLING {
    take:
    ch_assembly
    ch_reads
    ch_sample_id

    main:

    ch_versions = Channel.empty()

    if ( params.use_masking ) {
        // VARIANT CALLING
        bwa_index(ch_assembly)
        minimap2_index(ch_assembly)

        // create input map channels for bwa on assembly
        ch_reads
            .join(bwa_index.out.index)
            .multiMap { id, reads, index -> 
                reads: tuple(id, reads)
                index: index
            }
            .set{ ch_bwa_mem_assembly_map }
        bwa_mem_assembly(ch_bwa_mem_assembly_map.reads, ch_bwa_mem_assembly_map.index)

        // create input map channels for minimap2 on assembly
        ch_reads
            .join(minimap2_index.out.index)
            .multiMap { id, reads, index -> 
                reads: tuple(id, reads)
                index: index
            }
            .set{ ch_minimap2_align_assembly_map }
        minimap2_align_assembly(ch_minimap2_align_assembly_map.reads, ch_minimap2_align_assembly_map.index)
        samtools_sort_assembly(minimap2_align_assembly.out.sam)

        bwa_mem_assembly.out.bam.mix(samtools_sort_assembly.out.bam).set{ ch_bam }

        samtools_index_assembly(ch_bam)

        // construct freebayes input channels
        ch_bam
            .join(samtools_index_assembly.out.bai)
            .set{ ch_bam_bai }

        freebayes(ch_assembly.join(ch_bam_bai))
        ch_vcf = freebayes.out.vcf

        ch_versions = ch_versions.mix(bwa_index.out.versions)
        ch_versions = ch_versions.mix(bwa_mem_assembly.out.versions)
        ch_versions = ch_versions.mix(freebayes.out.versions)
        ch_versions = ch_versions.mix(minimap2_align_assembly.out.versions)
        ch_versions = ch_versions.mix(minimap2_index.out.versions)
        ch_versions = ch_versions.mix(samtools_index_assembly.out.versions)
        ch_versions = ch_versions.mix(samtools_sort_assembly.out.versions)

    } else {
      ch_sample_id.set{ ch_vcf }
    }

    emit:
    vcf         = ch_vcf                // channel: [ val(meta), path(vcf) ]
    versions    = ch_versions           // channel: [ versions.yml ]
}
