#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { post_align_qc         } from '../modules/local/prp/main.nf'
include { CALL_ASSEMBLY         } from '../subworkflows/assembly.nf'
include { CALL_POSTPROCESSING   } from '../subworkflows/postprocessing.nf'
include { CALL_PREPROCESSING    } from '../subworkflows/preprocessing.nf'
include { CALL_PROFILING        } from '../subworkflows/profiling.nf'
include { CALL_QUALITY_CONTROL  } from '../subworkflows/quality_control.nf'
include { CALL_RELATEDNESS      } from '../subworkflows/relatedness.nf'

workflow CALL_MYCOBACTERIUM_TUBERCULOSIS {
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
    gambit_db               = params.gambit_db ? file(params.gambit_db, checkIfExists: true) : Channel.value([])
    kraken_db               = params.kraken_db ? file(params.kraken_db, checkIfExists: true) : Channel.value([])
    mlst_blast_db           = params.mlst_blast_db ? file(params.mlst_blast_db, checkIfExists: true) : Channel.value([])
    pointfinder_db          = params.pointfinder_db ? file(params.pointfinder_db, checkIfExists: true) : Channel.value([])
    pubmlst_db              = params.pubmlst_db ? file(params.pubmlst_db, checkIfExists: true) : Channel.value([])
    resfinder_db            = params.resfinder_db ? file(params.resfinder_db, checkIfExists: true) : Channel.value([])
    serotypefinder_db       = params.serotypefinder_db ? file(params.serotypefinder_db, checkIfExists: true) : Channel.value([])
    shigapass_db            = params.shigapass_db ? file(params.shigapass_db, checkIfExists: true) : Channel.value([])
    tb_grading_rules_bed    = params.tb_grading_rules_bed ? file(params.tb_grading_rules_bed, checkIfExists: true) : Channel.value([])
    tbdb_bed                = params.tbdb_bed ? file(params.tbdb_bed, checkIfExists: true) : Channel.value([])
    tbdb_bed_idx            = params.tbdb_bed_idx ? file(params.tbdb_bed_idx, checkIfExists: true) : Channel.value([])
    training_file           = params.training_file ? file(params.training_file, checkIfExists: true) : Channel.value([])
    virulencefinder_db      = params.virulencefinder_db ? file(params.virulencefinder_db, checkIfExists: true) : Channel.value([])

    // schemas and values
    assay                   = params.assay ? params.assay : Channel.value([])
    hostile_dir             = params.hostile_dir ? file(params.hostile_dir, checkIfExists: true) : Channel.value([])
    hostile_idx             = params.hostile_idx ? params.hostile_idx : Channel.value([])
    mlst_scheme             = params.mlst_scheme ? params.mlst_scheme : Channel.value([])
    platform                = params.platform ? params.platform : Channel.value([])
    release_life_cycle      = params.release_life_cycle ? params.release_life_cycle : Channel.value([])
    species                 = params.species ? params.species : Channel.value([])
    species_dir             = params.species_dir ? params.species_dir : Channel.value([])
    target_sample_size      = params.target_sample_size ? params.target_sample_size : Channel.value([])

    main:
    ch_versions = Channel.empty()

    CALL_PREPROCESSING (
        assay,
        hostile_dir,
        hostile_idx,
        input_samples,
        platform,
        release_life_cycle
    )

    CALL_ASSEMBLY (
        CALL_PREPROCESSING.out.reads
    )

    CALL_QUALITY_CONTROL (
        core_loci_bed,
        gambit_db,
        kraken_db,
        reference_genome,
        reference_genome_dir,
        reference_genome_idx,
        CALL_ASSEMBLY.out.assembly,
        CALL_PREPROCESSING.out.empty,
        CALL_PREPROCESSING.out.reads
    )

    CALL_RELATEDNESS (
        CALL_ASSEMBLY.out.assembly,
        CALL_PREPROCESSING.out.reads
    )

    CALL_PROFILING (
        core_loci_bed,
        reference_genome,
        tbdb_bed,
        tbdb_bed_idx,
        CALL_ASSEMBLY.out.assembly,
        CALL_PREPROCESSING.out.reads
    )

    post_align_qc(CALL_PROFILING.out.bam, reference_genome, core_loci_bed)

    CALL_PROFILING.out.bam
        .join(CALL_PROFILING.out.bai)
        .join(CALL_QUALITY_CONTROL.out.gambitcore)
        .join(CALL_QUALITY_CONTROL.out.kraken)
        .join(post_align_qc.out.json)
        .join(CALL_QUALITY_CONTROL.out.quast)
        .set{ ch_qc_combined_output }

    CALL_PREPROCESSING.out.empty
        .map{ sample_id, empty -> [ sample_id, empty, empty, empty, empty, empty, empty, empty ] }
        .set{ ch_screening_combined_output }
    
    CALL_PREPROCESSING.out.empty
        .map{ sample_id, empty -> [ sample_id, empty, empty, empty, empty, empty, empty, empty, empty ] }
        .set{ ch_typing_combined_output }

    CALL_POSTPROCESSING (
        reference_genome,
        reference_genome_idx,
        reference_genome_gff,
        species_dir,
        tb_grading_rules_bed,
        tbdb_bed,
        CALL_PREPROCESSING.out.combined_output,
        CALL_PROFILING.out.combined_output,
        ch_qc_combined_output,
        CALL_RELATEDNESS.out.combined_output,
        ch_screening_combined_output,
        CALL_PREPROCESSING.out.seqrun_meta,
        ch_typing_combined_output,
        CALL_PROFILING.out.vcf
    )

    ch_versions = ch_versions.mix(CALL_ASSEMBLY.out.versions)
    ch_versions = ch_versions.mix(CALL_PREPROCESSING.out.versions)
    ch_versions = ch_versions.mix(CALL_PROFILING.out.versions)
    ch_versions = ch_versions.mix(CALL_POSTPROCESSING.out.versions)
    ch_versions = ch_versions.mix(CALL_QUALITY_CONTROL.out.versions)
    ch_versions = ch_versions.mix(CALL_RELATEDNESS.out.versions)

    emit: 
    pipeline_result = CALL_POSTPROCESSING.out.pipeline_result   // channel: [ path(json) ]
    cdm             = CALL_POSTPROCESSING.out.cdm               // channel: [ path(txt) ]
    yaml            = CALL_POSTPROCESSING.out.yaml              // channel: [ path(yaml) ]
    versions        = ch_versions                               // channel: [ versions.yml ]
}
