#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { get_meta                                  } from '../methods/get_sample_data.nf'
include { get_reads                                 } from '../methods/get_sample_data.nf'
include { amrfinderplus                             } from '../modules/nf-core/amrfinderplus/main.nf'
include { bracken                                   } from '../modules/nf-core/bracken/main.nf'
include { bwa_index                                 } from '../modules/nf-core/bwa/main.nf'
include { bwa_mem as bwa_mem_assembly               } from '../modules/nf-core/bwa/main.nf'
include { chewbbaca_allelecall                      } from '../modules/nf-core/chewbbaca/main.nf'
include { chewbbaca_create_batch_list               } from '../modules/local/chewbbaca/batch/main.nf'
include { chewbbaca_split_results                   } from '../modules/local/chewbbaca/split/main.nf'
include { create_analysis_result                    } from '../modules/local/prp/main.nf'
include { create_cdm_input                          } from '../modules/local/prp/main.nf'
include { create_prp_yaml                           } from '../modules/local/yaml/prp/main.nf'
include { emmtyper                                  } from '../modules/nf-core/emmtyper/main.nf'
include { export_to_cdm                             } from '../modules/local/cmd/main.nf'
include { freebayes                                 } from '../modules/nf-core/freebayes/main.nf'
include { kraken                                    } from '../modules/nf-core/kraken/main.nf'
include { mask_polymorph_assembly                   } from '../modules/local/mask/main.nf'
include { minimap2_align as minimap2_align_assembly } from '../modules/nf-core/minimap2/main.nf'       
include { minimap2_index                            } from '../modules/nf-core/minimap2/main.nf'       
include { mlst                                      } from '../modules/nf-core/mlst/main.nf'
include { resfinder                                 } from '../modules/nf-core/resfinder/main.nf'
include { samtools_index as samtools_index_assembly } from '../modules/nf-core/samtools/main.nf'
include { samtools_sort as samtools_sort_assembly   } from '../modules/nf-core/samtools/main.nf'
include { serotypefinder                            } from '../modules/nf-core/serotypefinder/main.nf'
include { shigapass                                 } from '../modules/nf-core/shigapass/main.nf'
include { virulencefinder                           } from '../modules/nf-core/virulencefinder/main.nf'
include { CALL_BACTERIAL_BASE                       } from '../workflows/bacterial_base.nf'

workflow CALL_STREPTOCOCCUS {
    // set input data
    input_samples       = file(params.csv, checkIfExists: true)

    // load references
    reference_genome     = params.reference_genome ? file(params.reference_genome, checkIfExists: true) : Channel.value([])
    reference_genome_dir = params.reference_genome ? file(reference_genome.getParent(), checkIfExists: true) : Channel.value([])
    reference_genome_gff = params.reference_genome_gff ? file(params.reference_genome_gff, checkIfExists: true) : Channel.value([])
    reference_genome_idx = params.reference_genome_idx ? file(params.reference_genome_idx, checkIfExists: true) : Channel.value([])

    // databases
    amrfinder_db        = params.amrfinder_db ? file(params.amrfinder_db, checkIfExists: true) : Channel.value([])
    chewbbaca_db        = params.chewbbaca_db ? file(params.chewbbaca_db, checkIfExists: true) : Channel.value([])
    core_loci_bed       = params.core_loci_bed ? file(params.core_loci_bed, checkIfExists: true) : Channel.value([])
    kraken_db           = params.kraken_db ? file(params.kraken_db, checkIfExists: true) : Channel.value([])
    mlst_blast_db       = params.mlst_blast_db ? file(params.mlst_blast_db, checkIfExists: true) : Channel.value([])
    pointfinder_db      = params.pointfinder_db ? file(params.pointfinder_db, checkIfExists: true) : Channel.value([])
    pubmlst_db          = params.pubmlst_db ? file(params.pubmlst_db, checkIfExists: true) : Channel.value([])
    resfinder_db        = params.resfinder_db ? file(params.resfinder_db, checkIfExists: true) : Channel.value([])
    serotypefinder_db   = params.serotypefinder_db ? file(params.serotypefinder_db, checkIfExists: true) : Channel.value([])
    shigapass_db        = params.shigapass_db ? file(params.shigapass_db, checkIfExists: true) : Channel.value([])
    training_file       = params.training_file ? file(params.training_file, checkIfExists: true) : Channel.value([])
    virulencefinder_db  = params.virulencefinder_db ? file(params.virulencefinder_db, checkIfExists: true) : Channel.value([])

    // schemas and values
    mlst_scheme         = params.mlst_scheme ? params.mlst_scheme : Channel.value([])
    species             = params.species ? params.species : Channel.value([])
    species_dir         = params.species_dir ? params.species_dir : Channel.value([])
    target_sample_size  = params.target_sample_size ? params.target_sample_size : Channel.value([])

    main:
        ch_versions = Channel.empty()

        CALL_BACTERIAL_BASE( core_loci_bed, reference_genome, reference_genome_dir, reference_genome_idx, input_samples, kraken_db, target_sample_size )
        
        CALL_BACTERIAL_BASE.out.assembly.set{ch_assembly}
        CALL_BACTERIAL_BASE.out.bam.set{ch_ref_bam}
        CALL_BACTERIAL_BASE.out.bai.set{ch_ref_bai}
        CALL_BACTERIAL_BASE.out.empty.set{ch_empty}
        CALL_BACTERIAL_BASE.out.id_meta.set{ch_id_meta}
        CALL_BACTERIAL_BASE.out.kraken.set{ch_kraken}
        CALL_BACTERIAL_BASE.out.metadata.set{ch_metadata}
        CALL_BACTERIAL_BASE.out.quast.set{ch_quast}
        CALL_BACTERIAL_BASE.out.qc.set{ch_qc}
        CALL_BACTERIAL_BASE.out.reads.set{ch_reads}
        CALL_BACTERIAL_BASE.out.reads_w_meta.set{ch_input_meta}
        CALL_BACTERIAL_BASE.out.seqplat_meta.set{ch_seqplat_meta}
        CALL_BACTERIAL_BASE.out.seqrun_meta.set{ch_seqrun_meta}
        CALL_BACTERIAL_BASE.out.ska_build.set{ch_ska}
        CALL_BACTERIAL_BASE.out.sourmash.set{ch_sourmash}

        bwa_index(ch_assembly.join(ch_seqplat_meta))
        minimap2_index(ch_assembly.join(ch_seqplat_meta))

        // create input map channels for bwa on assembly
        ch_input_meta
            .join(bwa_index.out.index)
            .multiMap { id, reads, platform, index -> 
                reads_w_meta: tuple(id, reads, platform)
                index: index
            }
            .set{ ch_bwa_mem_assembly_map }
        bwa_mem_assembly(ch_bwa_mem_assembly_map.reads_w_meta, ch_bwa_mem_assembly_map.index)

        // create input map channels for minimap2 on assembly
        ch_input_meta
            .join(minimap2_index.out.index)
            .multiMap { id, reads, platform, index -> 
                reads_w_meta: tuple(id, reads, platform)
                index: index
            }
            .set{ ch_minimap2_align_assembly_map }
        minimap2_align_assembly(ch_minimap2_align_assembly_map.reads_w_meta, ch_minimap2_align_assembly_map.index)
        samtools_sort_assembly(minimap2_align_assembly.out.sam)

        bwa_mem_assembly.out.bam.mix(samtools_sort_assembly.out.bam).set{ ch_bam }

        samtools_index_assembly(ch_bam)

        // construct freebayes input channels
        ch_bam
            .join(samtools_index_assembly.out.bai)
            .set{ ch_bam_bai }

        // VARIANT CALLING
        freebayes(ch_assembly, ch_bam_bai)

        mask_polymorph_assembly(ch_assembly.join(freebayes.out.vcf))

        // TYPING
        mlst(ch_assembly, mlst_scheme, pubmlst_db, mlst_blast_db)

        mask_polymorph_assembly.out.fasta
            .multiMap { sample_id, filePath -> 
                sample_id: sample_id
                filePath: filePath
            }
            .set{ maskedAssemblyMap }

        chewbbaca_create_batch_list(maskedAssemblyMap.filePath.collect())
        chewbbaca_allelecall(chewbbaca_create_batch_list.out.list, chewbbaca_db, training_file)
        chewbbaca_split_results(maskedAssemblyMap.sample_id.collect(), chewbbaca_allelecall.out.calls)
        emmtyper(ch_assembly)
        serotypefinder(ch_assembly, params.use_serotype_dbs, serotypefinder_db)
        shigapass(ch_assembly, shigapass_db)

        // SCREENING
        // antimicrobial detection (amrfinderplus)
        amrfinderplus(ch_assembly, species, amrfinder_db)

        // resistance & virulence prediction
        resfinder(ch_input_meta, species, resfinder_db, pointfinder_db)
        virulencefinder(ch_reads, params.use_virulence_dbs, virulencefinder_db)

        ch_reads.map{ sample_id, reads -> [ sample_id, [] ] }.set{ ch_empty }
        
        ch_id_meta
            .join(ch_quast)
            .join(ch_empty)
            .join(mlst.out.json)
            .join(chewbbaca_split_results.out.output)
            .join(amrfinderplus.out.output)
            .join(resfinder.out.json)
            .join(resfinder.out.meta)
            .join(ch_empty)
            .join(ch_empty)
            .join(virulencefinder.out.json)
            .join(virulencefinder.out.meta)
            .join(ch_empty)
            .join(emmtyper.out.tsv)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_metadata)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_kraken)
            .join(ch_sourmash)
            .join(ch_ska)
            .set{ combined_output }

        create_prp_yaml(combined_output, reference_genome, reference_genome_idx, reference_genome_gff)

        create_analysis_result(create_prp_yaml.out.yaml)

        ch_quast
            .join(ch_empty)
            .join(chewbbaca_split_results.out.output)
            .set{ ch_cdm_input }

        create_cdm_input(ch_cdm_input)

        export_to_cdm(create_cdm_input.out.json.join(ch_seqrun_meta), species_dir)

        ch_versions = ch_versions.mix(CALL_BACTERIAL_BASE.out.versions)
        ch_versions = ch_versions.mix(amrfinderplus.out.versions)
        ch_versions = ch_versions.mix(bracken.out.versions)
        ch_versions = ch_versions.mix(bwa_index.out.versions)
        ch_versions = ch_versions.mix(bwa_mem_assembly.out.versions)
        ch_versions = ch_versions.mix(chewbbaca_allelecall.out.versions)
        ch_versions = ch_versions.mix(create_analysis_result.out.versions)
        ch_versions = ch_versions.mix(emmtyper.out.versions)
        ch_versions = ch_versions.mix(freebayes.out.versions)
        ch_versions = ch_versions.mix(kraken.out.versions)
        ch_versions = ch_versions.mix(mlst.out.versions)
        ch_versions = ch_versions.mix(resfinder.out.versions)
        ch_versions = ch_versions.mix(samtools_index_assembly.out.versions)
        ch_versions = ch_versions.mix(serotypefinder.out.versions)
        ch_versions = ch_versions.mix(shigapass.out.versions)
        ch_versions = ch_versions.mix(virulencefinder.out.versions)

    emit: 
        pipeline_result = create_analysis_result.out.json
        cdm             = export_to_cdm.out.cdm
        yaml            = create_yaml.out.yaml
        versions        = ch_versions
}
