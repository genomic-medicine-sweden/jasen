#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { amrfinderplus                             } from '../nextflow-modules/modules/amrfinderplus/main.nf'
include { bracken                                   } from '../nextflow-modules/modules/bracken/main.nf'
include { bwa_index                                 } from '../nextflow-modules/modules/bwa/main.nf'
include { bwa_mem as bwa_mem_assembly               } from '../nextflow-modules/modules/bwa/main.nf'
include { chewbbaca_allelecall                      } from '../nextflow-modules/modules/chewbbaca/main.nf'
include { chewbbaca_create_batch_list               } from '../nextflow-modules/modules/chewbbaca/main.nf'
include { chewbbaca_split_results                   } from '../nextflow-modules/modules/chewbbaca/main.nf'
include { create_analysis_result                    } from '../nextflow-modules/modules/prp/main.nf'
include { create_cdm_input                          } from '../nextflow-modules/modules/prp/main.nf'
include { create_yaml                               } from '../nextflow-modules/modules/yaml/main.nf'
include { emmtyper                                  } from '../nextflow-modules/modules/emmtyper/main.nf'
include { export_to_cdm                             } from '../nextflow-modules/modules/cmd/main.nf'
include { freebayes                                 } from '../nextflow-modules/modules/freebayes/main.nf'
include { kraken                                    } from '../nextflow-modules/modules/kraken/main.nf'
include { mask_polymorph_assembly                   } from '../nextflow-modules/modules/mask/main.nf'
include { minimap2_align as minimap2_align_assembly } from '../nextflow-modules/modules/minimap2/main.nf'       
include { minimap2_index                            } from '../nextflow-modules/modules/minimap2/main.nf'      
include { mlst                                      } from '../nextflow-modules/modules/mlst/main.nf'
include { resfinder                                 } from '../nextflow-modules/modules/resfinder/main.nf'
include { samtools_index as samtools_index_assembly } from '../nextflow-modules/modules/samtools/main.nf'
include { samtools_sort as samtools_sort_assembly   } from '../nextflow-modules/modules/samtools/main.nf'
include { virulencefinder                           } from '../nextflow-modules/modules/virulencefinder/main.nf'
include { CALL_BACTERIAL_BASE                       } from '../workflows/bacterial_base.nf'

workflow CALL_STREPTOCOCCUS_PYOGENES {
    // set input data
    inputSamples = file(params.csv, checkIfExists: true)

    // load references
    referenceGenome = params.referenceGenome ? file(params.referenceGenome, checkIfExists: true) : Channel.value([])
    referenceGenomeDir = params.referenceGenome ? file(referenceGenome.getParent(), checkIfExists: true) : Channel.value([])
    referenceGenomeGff = params.referenceGenomeGff ? file(params.referenceGenomeGff, checkIfExists: true) : Channel.value([])
    referenceGenomeFai = params.referenceGenomeFai ? file(params.referenceGenomeFai, checkIfExists: true) : Channel.value([])
    referenceGenomeMmi = params.referenceGenomeMmi ? file(params.referenceGenomeMmi, checkIfExists: true) : Channel.value([])
    // databases
    amrfinderDb = file(params.amrfinderDb, checkIfExists: true)
    chewbbacaDb = file(params.chewbbacaDb, checkIfExists: true)
    coreLociBed = params.coreLociBed ? file(params.coreLociBed, checkIfExists: true) : Channel.value([])
    krakenDb = params.krakenDb ? file(params.krakenDb, checkIfExists: true) : Channel.value([])
    mlstBlastDb = params.mlstBlastDb ? file(params.mlstBlastDb, checkIfExists: true) : Channel.value([])
    pointfinderDb = file(params.pointfinderDb, checkIfExists: true)
    pubMlstDb = params.pubMlstDb ? file(params.pubMlstDb, checkIfExists: true) : Channel.value([])
    resfinderDb = file(params.resfinderDb, checkIfExists: true)
    trainingFile = params.trainingFile ? file(params.trainingFile, checkIfExists: true) : Channel.value([])
    virulencefinderDb = file(params.virulencefinderDb, checkIfExists: true)
    // schemas and values
    mlstScheme = params.mlstScheme ? params.mlstScheme : Channel.value([])
    species = params.species ? params.species : Channel.value([])
    speciesDir = params.speciesDir ? params.speciesDir : Channel.value([])
    targetSampleSize = params.targetSampleSize ? params.targetSampleSize : Channel.value([])

    main:
        ch_versions = Channel.empty()

        CALL_BACTERIAL_BASE( coreLociBed, referenceGenome, referenceGenomeDir, referenceGenomeMmi, inputSamples, targetSampleSize )
        
        CALL_BACTERIAL_BASE.out.assembly.set{ch_assembly}
        CALL_BACTERIAL_BASE.out.reads.set{ch_reads}
        CALL_BACTERIAL_BASE.out.bam.set{ch_ref_bam}
        CALL_BACTERIAL_BASE.out.bai.set{ch_ref_bai}
        CALL_BACTERIAL_BASE.out.quast.set{ch_quast}
        CALL_BACTERIAL_BASE.out.qc.set{ch_qc}
        CALL_BACTERIAL_BASE.out.metadata.set{ch_metadata}
        CALL_BACTERIAL_BASE.out.seqplat_meta.set{ch_seqplat_meta}
        CALL_BACTERIAL_BASE.out.seqrun_meta.set{ch_seqrun_meta}
        CALL_BACTERIAL_BASE.out.reads_w_meta.set{ch_input_meta}
        CALL_BACTERIAL_BASE.out.sourmash.set{ch_sourmash}
        CALL_BACTERIAL_BASE.out.ska_build.set{ch_ska}

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
        mlst(ch_assembly, params.mlstScheme, pubMlstDb, mlstBlastDb)

        mask_polymorph_assembly.out.fasta
            .multiMap { sampleID, filePath -> 
                sampleID: sampleID
                filePath: filePath
            }
            .set{ maskedAssemblyMap }

        chewbbaca_create_batch_list(maskedAssemblyMap.filePath.collect())
        chewbbaca_allelecall(chewbbaca_create_batch_list.out.list, chewbbacaDb, trainingFile)
        chewbbaca_split_results(maskedAssemblyMap.sampleID.collect(), chewbbaca_allelecall.out.calls)
        emmtyper(ch_assembly)

        // SCREENING
        // antimicrobial detection (amrfinderplus)
        amrfinderplus(ch_assembly, params.species, amrfinderDb)

        // resistance & virulence prediction
        resfinder(ch_input_meta, params.species, resfinderDb, pointfinderDb)
        virulencefinder(ch_reads, params.useVirulenceDbs, virulencefinderDb)

        ch_reads.map { sampleID, reads -> [ sampleID, [] ] }.set{ ch_empty }

        ch_quast
            .join(ch_qc)
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
            .join(ch_ref_bam)
            .join(ch_ref_bai)
            .join(ch_metadata)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .set{ combinedOutput }

        if ( params.useKraken ) {
            krakenDb = file(params.krakenDb, checkIfExists: true)
            kraken(ch_reads, krakenDb)
            bracken(kraken.out.report, krakenDb).output
            combinedOutput.join(bracken.out.output).set{ combinedOutput }
            create_analysis_result(combinedOutput, referenceGenome, referenceGenomeFai, referenceGenomeMmi, referenceGenomeGff)
            ch_versions = ch_versions.mix(kraken.out.versions)
            ch_versions = ch_versions.mix(bracken.out.versions)
        } else {
            combinedOutput.join(ch_empty).set{ combinedOutput }
            create_analysis_result(combinedOutput, referenceGenome, referenceGenomeFai, referenceGenomeMmi, referenceGenomeGff)
        }

        create_yaml(create_analysis_result.out.json.join(ch_sourmash).join(ch_ska), speciesDir)

        ch_quast
            .join(ch_qc)
            .join(chewbbaca_split_results.out.output)
            .set{ cdmInput }

        create_cdm_input(cdmInput)

        export_to_cdm(create_cdm_input.out.json.join(ch_seqrun_meta), params.speciesDir)

        ch_versions = ch_versions.mix(CALL_BACTERIAL_BASE.out.versions)
        ch_versions = ch_versions.mix(amrfinderplus.out.versions)
        ch_versions = ch_versions.mix(bwa_index.out.versions)
        ch_versions = ch_versions.mix(bwa_mem_assembly.out.versions)
        ch_versions = ch_versions.mix(chewbbaca_allelecall.out.versions)
        ch_versions = ch_versions.mix(create_analysis_result.out.versions)
        ch_versions = ch_versions.mix(freebayes.out.versions)
        ch_versions = ch_versions.mix(mlst.out.versions)
        ch_versions = ch_versions.mix(resfinder.out.versions)
        ch_versions = ch_versions.mix(samtools_index_assembly.out.versions)
        ch_versions = ch_versions.mix(virulencefinder.out.versions)

    emit: 
        pipeline_result = create_analysis_result.out.json
        cdm             = export_to_cdm.out.cdm
        yaml            = create_yaml.out.yaml
        versions        = ch_versions
}
