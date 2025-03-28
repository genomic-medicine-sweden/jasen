#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { get_meta                                  } from '../methods/get_sample_data.nf'
include { get_reads                                 } from '../methods/get_sample_data.nf'
include { annotate_delly                            } from '../modules/local/prp/main.nf'
include { bracken                                   } from '../modules/nf-core/bracken/main.nf'
include { create_analysis_result                    } from '../modules/local/prp/main.nf'
include { add_igv_track as add_tbdb_bed_track       } from '../modules/local/prp/main.nf'
include { add_igv_track as add_grading_bed_track    } from '../modules/local/prp/main.nf'
include { create_cdm_input                          } from '../modules/local/prp/main.nf'
include { create_prp_yaml                           } from '../modules/local/yaml/prp/main.nf'
include { export_to_cdm                             } from '../modules/local/cmd/main.nf'
include { kraken                                    } from '../modules/nf-core/kraken/main.nf'
include { mykrobe                                   } from '../modules/nf-core/mykrobe/main.nf'
include { post_align_qc                             } from '../modules/local/prp/main.nf'
include { snippy                                    } from '../modules/nf-core/snippy/main.nf'
include { tbprofiler as tbprofiler_mergedb          } from '../modules/nf-core/tbprofiler/main.nf'
include { CALL_BACTERIAL_BASE                       } from '../workflows/bacterial_base.nf'

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

        CALL_BACTERIAL_BASE( core_loci_bed, reference_genome, reference_genome_dir, input_samples, kraken_db, target_sample_size )

        CALL_BACTERIAL_BASE.out.assembly.set{ch_assembly}
        CALL_BACTERIAL_BASE.out.empty.set{ch_empty}
        CALL_BACTERIAL_BASE.out.id_meta.set{ch_id_meta}
        CALL_BACTERIAL_BASE.out.kraken.set{ch_kraken}
        CALL_BACTERIAL_BASE.out.metadata.set{ch_metadata}
        CALL_BACTERIAL_BASE.out.quast.set{ch_quast}
        CALL_BACTERIAL_BASE.out.reads.set{ch_reads}
        CALL_BACTERIAL_BASE.out.reads_w_meta.set{ch_input_meta}
        CALL_BACTERIAL_BASE.out.seqplat_meta.set{ch_seqplat_meta}
        CALL_BACTERIAL_BASE.out.seqrun_meta.set{ch_seqrun_meta}
        CALL_BACTERIAL_BASE.out.ska_build.set{ch_ska}
        CALL_BACTERIAL_BASE.out.sourmash.set{ch_sourmash}

        mykrobe(ch_reads)

        snippy(ch_reads, reference_genome)

        tbprofiler_mergedb(ch_reads)

        annotate_delly(tbprofiler_mergedb.out.vcf, tbdb_bed, tbdb_bed_idx)

        post_align_qc(tbprofiler_mergedb.out.bam, params.reference_genome, core_loci_bed)
        post_align_qc.out.qc.set{ch_qc}

        ch_id_meta
            .join(ch_quast)
            .join(ch_qc)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(tbprofiler_mergedb.out.bam)
            .join(tbprofiler_mergedb.out.bai)
            .join(ch_metadata)
            .join(annotate_delly.out.vcf)
            .join(mykrobe.out.csv)
            .join(tbprofiler_mergedb.out.json)
            .join(ch_kraken)
            .join(ch_sourmash)
            .join(ch_ska)
            .set{ combined_output }

        create_prp_yaml(combined_output, reference_genome, reference_genome_idx, reference_genome_gff)

        create_analysis_result(create_prp_yaml.out.yaml)

        ch_quast
            .join(ch_qc)
            .join(ch_empty)
            .set{ ch_cdm_input }

        create_cdm_input(ch_cdm_input)

        export_to_cdm(create_cdm_input.out.json.join(ch_seqrun_meta), params.species_dir)

        ch_versions = ch_versions.mix(CALL_BACTERIAL_BASE.out.versions)
        ch_versions = ch_versions.mix(create_analysis_result.out.versions)
        ch_versions = ch_versions.mix(mykrobe.out.versions)
        ch_versions = ch_versions.mix(snippy.out.versions)
        ch_versions = ch_versions.mix(tbprofiler_mergedb.out.versions)

    emit: 
        pipeline_result = add_grading_bed_track.out.json
        cdm             = export_to_cdm.out.cdm
        yaml            = create_prp_yaml.out.yaml
        versions        = ch_versions
}
