include { bwa_mem as bwa_mem_ref                } from '../nextflow-modules/modules/bwa/main.nf'
include { samtools_index as samtools_index_ref  } from '../nextflow-modules/modules/samtools/main.nf'
include { post_align_qc                         } from '../nextflow-modules/modules/qc/main.nf'

workflow CALL_QUALITY_CONTROL {
    take:
        coreLociBed
        genomeReferenceDir
        ch_reads            // channel: [ val(meta), val(reads) ]

    main:
        ch_versions = Channel.empty()

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

        ch_versions = ch_versions.mix(bwa_mem_ref.out.versions)
        ch_versions = ch_versions.mix(samtools_index_ref.out.versions)

    emit:
        qc          = post_align_qc.out.qc  // channel: [ val(meta), path(fasta)]
        versions    = ch_versions           // channel: [ versions.yml ]
}