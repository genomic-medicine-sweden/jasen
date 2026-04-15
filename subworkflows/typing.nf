#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { bwa_index                                 } from '../modules/nf-core/bwa/main.nf'
include { bwa_mem as bwa_mem_assembly               } from '../modules/nf-core/bwa/main.nf'
include { chewbbaca_allelecall                      } from '../modules/nf-core/chewbbaca/main.nf'
include { chewbbaca_split_results                   } from '../modules/local/chewbbaca/split/main.nf'
include { emmtyper                                  } from '../modules/nf-core/emmtyper/main.nf'
include { freebayes                                 } from '../modules/nf-core/freebayes/main.nf'
include { minimap2_align as minimap2_align_assembly } from '../modules/nf-core/minimap2/main.nf'
include { minimap2_index                            } from '../modules/nf-core/minimap2/main.nf'
include { mlst                                      } from '../modules/nf-core/mlst/main.nf'
include { mask_polymorph_assembly                   } from '../modules/local/mask/main.nf'
include { samtools_index as samtools_index_assembly } from '../modules/nf-core/samtools/main.nf'
include { samtools_sort as samtools_sort_assembly   } from '../modules/nf-core/samtools/main.nf'
include { sccmec                                    } from '../modules/nf-core/sccmec/main.nf'
include { serotypefinder                            } from '../modules/nf-core/serotypefinder/main.nf'
include { shigapass                                 } from '../modules/nf-core/shigapass/main.nf'
include { spatyper                                  } from '../modules/nf-core/spatyper/main.nf'

workflow CALL_TYPING {
    take:
    chewbbaca_db
    mlst_blast_db
    mlst_scheme
    pubmlst_db
    serotypefinder_db
    shigapass_db
    species
    training_file
    ch_assembly
    ch_reads
    ch_sample_id

    main:

    ch_versions = Channel.empty()

    // TYPING
    if ( !params.ci ) {
        mlst(ch_assembly, mlst_scheme, pubmlst_db, mlst_blast_db)
        mlst.out.json.set{ ch_mlst }
        ch_versions = ch_versions.mix(mlst.out.versions)
    } else {
        ch_sample_id.set{ ch_mlst }
    }

    if ( params.use_masking ) {
        // MASKING VARIANT CALLING (against assembly, for cgMLST)
        bwa_index(ch_assembly)
        minimap2_index(ch_assembly)

        ch_reads
            .join(bwa_index.out.index)
            .multiMap { id, reads, index ->
                reads: tuple(id, reads)
                index: index
            }
            .set{ ch_bwa_mem_assembly_map }
        bwa_mem_assembly(ch_bwa_mem_assembly_map.reads, ch_bwa_mem_assembly_map.index)

        ch_reads
            .join(minimap2_index.out.index)
            .multiMap { id, reads, index ->
                reads: tuple(id, reads)
                index: index
            }
            .set{ ch_minimap2_align_assembly_map }
        minimap2_align_assembly(ch_minimap2_align_assembly_map.reads, ch_minimap2_align_assembly_map.index)
        samtools_sort_assembly(minimap2_align_assembly.out.sam)

        bwa_mem_assembly.out.bam.mix(samtools_sort_assembly.out.bam).set{ ch_assembly_bam }

        samtools_index_assembly(ch_assembly_bam)

        ch_assembly_bam
            .join(samtools_index_assembly.out.bai)
            .set{ ch_assembly_bam_bai }

        freebayes(ch_assembly.join(ch_assembly_bam_bai)).vcf.set{ ch_polymorph_vcf }

        ch_versions = ch_versions.mix(bwa_index.out.versions)
        ch_versions = ch_versions.mix(bwa_mem_assembly.out.versions)
        ch_versions = ch_versions.mix(freebayes.out.versions)
        ch_versions = ch_versions.mix(minimap2_align_assembly.out.versions)
        ch_versions = ch_versions.mix(minimap2_index.out.versions)
        ch_versions = ch_versions.mix(samtools_index_assembly.out.versions)
        ch_versions = ch_versions.mix(samtools_sort_assembly.out.versions)

        mask_polymorph_assembly(ch_assembly.join(ch_polymorph_vcf))

        mask_polymorph_assembly.out.fasta
            .multiMap { sample_id, fasta ->
                sample_id: sample_id
                fasta: fasta
            }
            .set{ assembly_map }
    } else {
        assembly_map = ch_assembly.multiMap { sample_id, fasta ->
            sample_id: sample_id
            fasta: fasta
        }
    }

    assembly_map.fasta
        .collectFile(name: 'batch_input.list', newLine: true) { fasta ->
            fasta.toRealPath().toString()
        }
        .set { batch_list }

    chewbbaca_allelecall(batch_list, chewbbaca_db, training_file)

    chewbbaca_split_results(assembly_map.sample_id.collect(), chewbbaca_allelecall.out.calls)

    // Call species-specific typing
    // ecoli
    if (params.species == "escherichia coli") {
        serotypefinder(ch_assembly, params.use_serotype_dbs, serotypefinder_db)
        serotypefinder.out.json.set{ ch_serotypefinder }
        serotypefinder.out.meta.set{ ch_serotypefinder_meta }
        shigapass(ch_assembly, shigapass_db).csv.set{ ch_shigapass }
        ch_versions = ch_versions.mix(serotypefinder.out.versions)
        ch_versions = ch_versions.mix(shigapass.out.versions)
    } else {
        ch_sample_id.set{ ch_serotypefinder }
        ch_sample_id.set{ ch_serotypefinder_meta }
        ch_sample_id.set{ ch_shigapass }
    }

    // saureus
    if (params.species == "staphylococcus aureus") {
        sccmec(ch_assembly).tsv.set{ ch_sccmec }
        spatyper(ch_assembly).tsv.set{ ch_spatyper }
        ch_versions = ch_versions.mix(sccmec.out.versions)
        ch_versions = ch_versions.mix(spatyper.out.versions)
    } else {
        ch_sample_id.set{ ch_sccmec }
        ch_sample_id.set{ ch_spatyper }
    }

    // streptococcus & spyogenes
    if (params.species in ["streptococcus", "streptococcus pyogenes"]) {
        emmtyper(ch_assembly).tsv.set{ ch_emmtyper }
        ch_versions = ch_versions.mix(emmtyper.out.versions)
    } else {
        ch_sample_id.set{ ch_emmtyper }
    }

    chewbbaca_split_results.out.tsv
        .join(ch_emmtyper)
        .join(ch_mlst)
        .join(ch_sccmec)
        .join(ch_serotypefinder)
        .join(ch_serotypefinder_meta)
        .join(ch_shigapass)
        .join(ch_spatyper)
        .set{ ch_combined_output }

    ch_versions = ch_versions.mix(chewbbaca_allelecall.out.versions)

    emit:
    chewbbaca       = chewbbaca_split_results.out.tsv     // channel: [ val(meta), path(tsv) ]
    combined_output = ch_combined_output                  // channel: [ val(meta), path(tsv), path(tsv), path(json), path(tsv), path(json), path(json), path(csv), path(tsv) ]
    emmtyper        = ch_emmtyper                         // channel: [ val(meta), path(tsv) ]
    mlst            = ch_mlst                             // channel: [ val(meta), path(json) ]
    sccmec          = ch_sccmec                           // channel: [ val(meta), path(tsv) ]
    serotypefinder  = ch_serotypefinder                   // channel: [ val(meta), path(json) ]
    serotypefinder  = ch_serotypefinder_meta              // channel: [ val(meta), path(json) ]
    shigapass       = ch_shigapass                        // channel: [ val(meta), path(csv) ]
    spatyper        = ch_spatyper                         // channel: [ val(meta), path(tsv) ]
    versions        = ch_versions                         // channel: [ versions.yml ]
}
