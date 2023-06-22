include { bwa_index                                 } from '../nextflow-modules/modules/bwa/main'
include { bwa_mem as bwa_mem_dedup                  } from '../nextflow-modules/modules/bwa/main'
include { freebayes                                 } from '../nextflow-modules/modules/freebayes/main'
include { samtools_index as samtools_index_assembly } from '../nextflow-modules/modules/samtools/main'

workflow CALL_VARIANT_CALLING {
    take:
        ch_assembly // channel: [ val(meta), val(fasta) ]
        ch_reads    // channel: [ val(meta), val(reads) ]

    main:
        ch_versions = Channel.empty()

        // mask polymorph regions
        bwa_index(ch_assembly)

        ch_reads
            .join(bwa_index.out.idx)
            .multiMap { id, reads, bai -> 
                reads: tuple(id, reads)
                bai: bai
            }
            .set { bwa_mem_dedup_ch }
        bwa_mem_dedup(bwa_mem_dedup_ch.reads, bwa_mem_dedup_ch.bai)
        samtools_index_assembly(bwa_mem_dedup.out.bam)

        // construct freebayes input channels
        ch_assembly
        .join(bwa_mem_dedup.out.bam)
        .join(samtools_index_assembly.out.bai)
        .multiMap { id, fasta, bam, bai -> 
            assembly: tuple(id, fasta)
            mapping: tuple(bam, bai)
        }
        .set { freebayes_ch }

        freebayes(freebayes_ch.assembly, freebayes_ch.mapping)

        ch_versions = ch_versions.mix(bwa_index.out.versions)
        ch_versions = ch_versions.mix(bwa_mem_dedup.out.versions)
        ch_versions = ch_versions.mix(freebayes.out.versions)
        ch_versions = ch_versions.mix(samtools_index_assembly.out.versions)

    emit:
        vcf         = freebayes.out.vcf // channel: [ val(meta), path(vcf)]
        versions    = ch_versions       // channel: [ versions.yml ]
}