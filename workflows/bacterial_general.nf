#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CALL_ASSEMBLY         } from '../subworkflows/assembly.nf'
include { CALL_POSTPROCESSING   } from '../subworkflows/postprocessing.nf'
include { CALL_PREPROCESSING    } from '../subworkflows/preprocessing.nf'
include { CALL_QUALITY_CONTROL  } from '../subworkflows/quality_control.nf'
include { CALL_RELATEDNESS      } from '../subworkflows/relatedness.nf'
include { CALL_SCREENING        } from '../subworkflows/screening.nf'
include { CALL_TYPING           } from '../subworkflows/typing.nf'
include { CALL_VARIANT_CALLING  } from '../subworkflows/variant_calling.nf'

workflow CALL_BACTERIAL_GENERAL {
    // set input data
    input_samples           = file(params.csv, checkIfExists: true)

    // load references 
    reference_genome        = params.reference_genome ? file(params.reference_genome, checkIfExists: true) : Channel.value([])
    reference_genome_dir    = params.reference_genome ? file(reference_genome.getParent(), checkIfExists: true) : Channel.value([])
    reference_genome_gff    = params.reference_genome_gff ? file(params.reference_genome_gff, checkIfExists: true) : Channel.value([])
    reference_genome_idx    = params.reference_genome_idx ? file(params.reference_genome_idx, checkIfExists: true) : Channel.value([])

    // databases
    amrfinder_db            = params.amrfinder_db ? file(params.amrfinder_db, checkIfExists: true) : Channel.value([])
    chewbbaca_db            = params.chewbbaca_db ? file(params.chewbbaca_db, checkIfExists: true) : Channel.value([])
    core_loci_bed           = params.core_loci_bed ? file(params.core_loci_bed, checkIfExists: true) : Channel.value([])
    kraken_db               = params.kraken_db ? file(params.kraken_db, checkIfExists: true) : Channel.value([])
    mlst_blast_db           = params.mlst_blast_db ? file(params.mlst_blast_db, checkIfExists: true) : Channel.value([])
    pointfinder_db          = params.pointfinder_db ? file(params.pointfinder_db, checkIfExists: true) : Channel.value([])
    pubmlst_db              = params.pubmlst_db ? file(params.pubmlst_db, checkIfExists: true) : Channel.value([])
    resfinder_db            = params.resfinder_db ? file(params.resfinder_db, checkIfExists: true) : Channel.value([])
    serotypefinder_db       = params.serotypefinder_db ? file(params.serotypefinder_db, checkIfExists: true) : Channel.value([])
    shigapass_db            = params.shigapass_db ? file(params.shigapass_db, checkIfExists: true) : Channel.value([])
    training_file           = params.training_file ? file(params.training_file, checkIfExists: true) : Channel.value([])
    virulencefinder_db      = params.virulencefinder_db ? file(params.virulencefinder_db, checkIfExists: true) : Channel.value([])

    // schemas and values
    hostile_dir             = params.hostile_dir ? file(params.hostile_dir, checkIfExists: true) : Channel.value([])
    hostile_idx             = params.hostile_idx ? params.hostile_idx : Channel.value([])
    mlst_scheme             = params.mlst_scheme ? params.mlst_scheme : Channel.value([])
    species                 = params.species ? params.species : Channel.value([])
    species_dir             = params.species_dir ? params.species_dir : Channel.value([])
    target_sample_size      = params.target_sample_size ? params.target_sample_size : Channel.value([])

    main:
    ch_versions = Channel.empty()

    CALL_PREPROCESSING (
        input_samples,
        hostile_dir,
        hostile_idx
    )

    CALL_ASSEMBLY (
        CALL_PREPROCESSING.out.reads_w_meta
    )

    CALL_QUALITY_CONTROL (
        core_loci_bed,
        kraken_db,
        reference_genome,
        reference_genome_dir,
        reference_genome_idx,
        CALL_ASSEMBLY.out.assembly,
        CALL_PREPROCESSING.out.empty,
        CALL_PREPROCESSING.out.reads,
        CALL_PREPROCESSING.out.reads_w_meta
    )

    CALL_RELATEDNESS (
        CALL_ASSEMBLY.out.assembly,
        CALL_PREPROCESSING.out.reads
    )

    CALL_VARIANT_CALLING (
        CALL_ASSEMBLY.out.assembly,
        CALL_PREPROCESSING.out.reads_w_meta,
        CALL_PREPROCESSING.out.reads,
        CALL_PREPROCESSING.out.seqplat_meta
    )

    CALL_TYPING (
        chewbbaca_db,
        mlst_blast_db,
        mlst_scheme,
        pubmlst_db,
        serotypefinder_db,
        shigapass_db,
        species,
        training_file,
        CALL_ASSEMBLY.out.assembly,
        CALL_PREPROCESSING.out.empty,
        CALL_VARIANT_CALLING.out.vcf
    )

    CALL_SCREENING (
        amrfinder_db,
        pointfinder_db,
        resfinder_db,
        virulencefinder_db,
        CALL_ASSEMBLY.out.assembly,
        CALL_PREPROCESSING.out.reads,
        CALL_PREPROCESSING.out.reads_w_meta
    )

    CALL_PREPROCESSING.out.empty
        .map{ sample_id, empty -> [ sample_id, empty, empty ] }
        .set{ ch_profiling_combined_output }

    CALL_POSTPROCESSING (
        reference_genome,
        reference_genome_idx,
        reference_genome_gff,
        species_dir,
        CALL_TYPING.out.chewbbaca,
        CALL_QUALITY_CONTROL.out.post_align_qc,
        CALL_PREPROCESSING.out.combined_output,
        ch_profiling_combined_output,
        CALL_QUALITY_CONTROL.out.combined_output,
        CALL_QUALITY_CONTROL.out.quast,
        CALL_RELATEDNESS.out.combined_output,
        CALL_SCREENING.out.combined_output,
        CALL_PREPROCESSING.out.seqrun_meta,
        CALL_TYPING.out.combined_output,
        CALL_VARIANT_CALLING.out.vcf
    )

    ch_versions = ch_versions.mix(CALL_ASSEMBLY.out.versions)
    ch_versions = ch_versions.mix(CALL_PREPROCESSING.out.versions)
    ch_versions = ch_versions.mix(CALL_POSTPROCESSING.out.versions)
    ch_versions = ch_versions.mix(CALL_QUALITY_CONTROL.out.versions)
    ch_versions = ch_versions.mix(CALL_RELATEDNESS.out.versions)
    ch_versions = ch_versions.mix(CALL_SCREENING.out.versions)
    ch_versions = ch_versions.mix(CALL_TYPING.out.versions)
    ch_versions = ch_versions.mix(CALL_VARIANT_CALLING.out.versions)

    emit:
    pipeline_result = CALL_POSTPROCESSING.out.pipeline_result   // channel: [ path(json) ]
    cdm             = CALL_POSTPROCESSING.out.cdm               // channel: [ path(txt) ]
    yaml            = CALL_POSTPROCESSING.out.yaml              // channel: [ path(yaml) ]
    versions        = ch_versions                               // channel: [ versions.yml ]
}
