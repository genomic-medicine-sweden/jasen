include { sourmash } from '../nextflow-modules/modules/sourmash/main.nf'

workflow CALL_RELATEDNESS {
    take:
        ch_assembly // channel: [ val(meta), val(fasta) ]

    main:
        ch_versions = Channel.empty()

        // sourmash
        sourmash(ch_assembly)
    
        ch_versions = ch_versions.mix(sourmash.out.versions)

    emit:
        sourmash    = sourmash.out.signature    // channel: [ val(meta), path(signature)]
        versions    = ch_versions               // channel: [ versions.yml ]
}