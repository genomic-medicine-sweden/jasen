include { quast             } from '../nextflow-modules/modules/quast/main'
include { skesa             } from '../nextflow-modules/modules/skesa/main'
include { spades_illumina   } from '../nextflow-modules/modules/spades/main'
include { spades_iontorrent } from '../nextflow-modules/modules/spades/main'

workflow CALL_ASSEMBLY {
    take:
        genomeReference
        ch_input_meta // channel: [ val(meta), val(input_meta) ]

    main:
        ch_versions = Channel.empty()

        // assembly
        skesa(ch_input_meta)
        spades_illumina(ch_input_meta)
        spades_iontorrent(ch_input_meta)

        Channel.empty().mix(skesa.out.fasta, spades_illumina.out.fasta, spades_iontorrent.out.fasta).set{ ch_assembly }

        // evaluate assembly quality 
        quast(ch_assembly, genomeReference)

        ch_versions = ch_versions.mix(skesa.out.versions)
        ch_versions = ch_versions.mix(spades_illumina.out.versions)
        ch_versions = ch_versions.mix(spades_iontorrent.out.versions)
        ch_versions = ch_versions.mix(quast.out.versions)

    emit:
        assembly    = ch_assembly   // channel: [ val(meta), path(fasta)]
        qc          = quast.out.qc  // channel: [ val(meta), path(qc)]
        versions    = ch_versions   // channel: [ versions.yml ]

}