#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { chewbbaca_allelecall          } from '../modules/nf-core/chewbbaca/main.nf'
include { chewbbaca_create_batch_list   } from '../modules/local/chewbbaca/batch/main.nf'
include { chewbbaca_split_results       } from '../modules/local/chewbbaca/split/main.nf'
include { emmtyper                      } from '../modules/nf-core/emmtyper/main.nf'
include { mlst                          } from '../modules/nf-core/mlst/main.nf'
include { mask_polymorph_assembly       } from '../modules/local/mask/main.nf'
include { sccmec                        } from '../modules/nf-core/sccmec/main.nf'
include { serotypefinder                } from '../modules/nf-core/serotypefinder/main.nf'
include { shigapass                     } from '../modules/nf-core/shigapass/main.nf'
include { spatyper                      } from '../modules/nf-core/spatyper/main.nf'

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
    ch_empty
    ch_vcf

    main:

    ch_versions = Channel.empty()
    ch_chewbbaca_input = Channel.empty()

    // TYPING
    if ( !params.ci ) {
        mlst(ch_assembly, mlst_scheme, pubmlst_db, mlst_blast_db)
        mlst.out.json.set{ ch_mlst }
        ch_versions = ch_versions.mix(mlst.out.versions)
    } else {
        ch_empty.set{ ch_mlst }
    }

    if ( params.use_masking ) {
        mask_polymorph_assembly(ch_assembly.join(ch_vcf))

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

    chewbbaca_create_batch_list(assembly_map.fasta.collect())
    chewbbaca_allelecall(chewbbaca_create_batch_list.out.list, chewbbaca_db, training_file)
    chewbbaca_split_results(assembly_map.sample_id.collect(), chewbbaca_allelecall.out.calls)

    // Call species-specific typing
    // ecoli
    if (params.species == "escherichia coli") {
        // serot
        serotypefinder(ch_assembly, params.use_serotype_dbs, serotypefinder_db)
        serotypefinder.out.json.set{ ch_serotypefinder }
        serotypefinder.out.meta.set{ ch_serotypefinder_meta }
        shigapass(ch_assembly, shigapass_db).csv.set{ ch_shigapass }
        ch_versions = ch_versions.mix(serotypefinder.out.versions)
        ch_versions = ch_versions.mix(shigapass.out.versions)
    } else {
        ch_empty.set{ ch_serotypefinder }
        ch_empty.set{ ch_serotypefinder_meta }
        ch_empty.set{ ch_shigapass }
    }

    // saureus
    if (params.species == "staphylococcus aureus") {
        sccmec(ch_assembly).tsv.set{ ch_sccmec }
        spatyper(ch_assembly).tsv.set{ ch_spatyper }
        ch_versions = ch_versions.mix(sccmec.out.versions)
        ch_versions = ch_versions.mix(spatyper.out.versions)
    } else {
        ch_empty.set{ ch_sccmec }
        ch_empty.set{ ch_spatyper }
    }

    // strep
    if (params.species in ["streptococcus", "streptococcus pyogenes"]) {
        emmtyper(ch_assembly).tsv.set{ ch_emmtyper }
        ch_versions = ch_versions.mix(emmtyper.out.versions)
    } else {
        ch_empty.set{ ch_emmtyper }
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
