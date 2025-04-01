#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { post_align_qc         } from '../modules/local/prp/main.nf'
include { CALL_ASSEMBLY         } from '../subworkflows/assembly.nf'
include { CALL_POSTPROCESSING   } from '../subworkflows/postprocessing.nf'
include { CALL_PREPROCESSING    } from '../subworkflows/preprocessing.nf'
include { CALL_PROFILING        } from '../subworkflows/profiling.nf'
include { CALL_QUALITY_CONTROL  } from '../subworkflows/quality_control.nf'

workflow CALL_MYCOBACTERIUM_TUBERCULOSIS {
    // set input data
    input_samples        = file(params.csv, checkIfExists: true)

    // load references 
    reference_genome     = file(params.reference_genome, checkIfExists: true)
    reference_genome_dir = file(reference_genome.getParent(), checkIfExists: true)
    reference_genome_idx = file(params.reference_genome_idx, checkIfExists: true)
    reference_genome_gff = file(params.reference_genome_gff, checkIfExists: true)

    // databases
    core_loci_bed        = file(params.core_loci_bed, checkIfExists: true)
    tbdb_bed             = file(params.tbdb_bed, checkIfExists: true)
    tbdb_bed_idx         = file(params.tbdb_bed_idx, checkIfExists: true)
    tb_grading_rules_bed = file(params.tb_grading_rules_bed, checkIfExists: true)

    // schemas and values
    target_sample_size   = params.target_sample_size ? params.target_sample_size : Channel.value([])

    main:
    ch_versions = Channel.empty()

    CALL_PREPROCESSING (
        input_samples
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

    CALL_PROFILING (
        core_loci_bed,
        reference_genome,
        tbdb_bed,
        tbdb_bed_idx,
        ch_assembly,
        ch_reads
    )

    post_align_qc(tbprofiler_mergedb.out.bam, reference_genome, core_loci_bed)

    CALL_PROFILING.out.bam
        .join(CALL_PROFILING.out.bai)
        .join(CALL_QUALITY_CONTROL.out.kraken)
        .join(post_align_qc.out.json)
        .join(CALL_QUALITY_CONTROL.out.quast)
        .set{ ch_qc_combined_output }

    CALL_PREPROCESSING.out.empty
        .map{ sample_id, empty -> [ sample_id, empty, empty, empty, empty, empty ] }
        .set{ ch_screening_combined_output }
    
    CALL_PREPROCESSING.out.empty
        .map{ sample_id, empty -> [ sample_id, empty, empty, empty, empty, empty, empty, empty, empty ] }
        .set{ ch_typing_combined_output }

    CALL_POSTPROCESSING (
        reference_genome
        reference_genome_idx,
        reference_genome_gff,
        species_dir,
        CALL_PREPROCESSING.out.empty,
        post_align_qc.out.json,
        CALL_PREPROCESSING.out.combined_output,
        CALL_PROFILING.out.combined_output,
        ch_qc_combined_output,
        CALL_QUALITY_CONTROL.out.quast,
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
    pipeline_result = add_grading_bed_track.out.json
    cdm             = export_to_cdm.out.cdm
    yaml            = create_prp_yaml.out.yaml
    versions        = ch_versions
}
