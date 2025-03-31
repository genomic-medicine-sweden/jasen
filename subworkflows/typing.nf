#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { chewbbaca_allelecall          } from '../modules/nf-core/chewbbaca/main.nf'
include { chewbbaca_create_batch_list   } from '../modules/local/chewbbaca/main.nf'
include { chewbbaca_split_results       } from '../modules/local/chewbbaca/main.nf'
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
    ch_versions

    main:
    // TYPING
    mlst(ch_assembly, mlst_scheme, pubmlst_db, mlst_blast_db)

    mask_polymorph_assembly(ch_assembly.join(ch_vcf))

    mask_polymorph_assembly.out.fasta
        .multiMap { sample_id, fasta -> 
            sample_id: sample_id
            fasta: fasta
        }
        .set{ masked_assembly_map }

    chewbbaca_create_batch_list(masked_assembly_map.fasta.collect())
    chewbbaca_allelecall(chewbbaca_create_batch_list.out.list, chewbbaca_db, training_file)
    chewbbaca_split_results(masked_assembly_map.sample_id.collect(), chewbbaca_allelecall.out.calls)

    // Call species-specific typing
    // ecoli
    if (params.species == "escherichia coli") {
        serotypefinder(ch_assembly, params.use_serotype_dbs, serotypefinder_db).set{ ch_serotypefinder }
        shigapass(ch_assembly, shigapass_db).set{ ch_shigapass }
        ch_versions = ch_versions.mix(serotypefinder.out.versions)
        ch_versions = ch_versions.mix(shigapass.out.versions)
    } else {
        ch_empty.set{ ch_serotypefinder }
        ch_empty.set{ ch_shigapass }
    }

    // saureus
    if (params.species == "staphylococcus aureus") {
        sccmec(ch_assembly).set{ ch_sccmec }
        spatyper(ch_assembly).set{ ch_spatyper }
        ch_versions = ch_versions.mix(sccmec.out.versions)
        ch_versions = ch_versions.mix(spatyper.out.versions)
    } else {
        ch_empty.set{ ch_sccmec }
        ch_empty.set{ ch_spatyper }
    }

    // strep
    if (params.species in ["streptococcus", "streptococcus pyogenes"]) {
        emmtyper(ch_assembly).set{ ch_emmtyper }
        ch_versions = ch_versions.mix(emmtyper.out.versions)
    } else {
        ch_empty.set{ ch_emmtyper }
    }

    ch_versions = ch_versions.mix(mlst.out.versions)
    ch_versions = ch_versions.mix(chewbbaca_allelecall.out.versions)

    emit:
    chewbbaca       = chewbbaca_split_results.out.tsv     // channel: [ val(meta), path(tsv)]
    mlst            = mlst.out.json                       // channel: [ val(meta), path(json)]
    sccmec          = ch_sccmec                           // channel: [ val(meta), path(tsv)]
    serotypefinder  = ch_serotypefinder                   // channel: [ val(meta), path(json)]
    shigapass       = ch_shigapass                        // channel: [ val(meta), path(csv)]
    spatyper        = ch_spatyper                         // channel: [ val(meta), path(tsv)]
    versions        = ch_versions                         // channel: [ versions.yml ]
}
