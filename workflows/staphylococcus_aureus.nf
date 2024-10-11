#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { get_meta                                  } from '../methods/get_meta.nf'
include { amrfinderplus                             } from '../nextflow-modules/modules/amrfinderplus/main.nf'
include { bracken                                   } from '../nextflow-modules/modules/bracken/main.nf'
include { bwa_index                                 } from '../nextflow-modules/modules/bwa/main.nf'
include { bwa_mem as bwa_mem_dedup                  } from '../nextflow-modules/modules/bwa/main.nf'
include { chewbbaca_allelecall                      } from '../nextflow-modules/modules/chewbbaca/main.nf'
include { chewbbaca_create_batch_list               } from '../nextflow-modules/modules/chewbbaca/main.nf'
include { chewbbaca_split_results                   } from '../nextflow-modules/modules/chewbbaca/main.nf'
include { copy_to_cron                              } from '../nextflow-modules/modules/cron/main.nf'
include { create_analysis_result                    } from '../nextflow-modules/modules/prp/main.nf'
include { create_cdm_input                          } from '../nextflow-modules/modules/prp/main.nf'
include { create_yaml                               } from '../nextflow-modules/modules/yaml/main.nf'
include { export_to_cdm                             } from '../nextflow-modules/modules/cmd/main.nf'
include { freebayes                                 } from '../nextflow-modules/modules/freebayes/main.nf'
include { kraken                                    } from '../nextflow-modules/modules/kraken/main.nf'
include { mask_polymorph_assembly                   } from '../nextflow-modules/modules/mask/main.nf'
include { mlst                                      } from '../nextflow-modules/modules/mlst/main.nf'
include { resfinder                                 } from '../nextflow-modules/modules/resfinder/main.nf'
include { samtools_index as samtools_index_assembly } from '../nextflow-modules/modules/samtools/main.nf'
include { serotypefinder                            } from '../nextflow-modules/modules/serotypefinder/main.nf'
include { virulencefinder                           } from '../nextflow-modules/modules/virulencefinder/main.nf'
include { CALL_BACTERIAL_BASE                       } from '../workflows/bacterial_base.nf'

workflow CALL_STAPHYLOCOCCUS_AUREUS {
    Channel.fromPath(params.csv).splitCsv(header:true)
        .map{ row -> get_meta(row) }
        .branch {
        iontorrent: it[2] == "iontorrent"
        illumina: it[2] == "illumina"
        nanopore: it[2] == "nanopore"
        }
        .set{ ch_meta }

    // load references 
    referenceGenome = file(params.referenceGenome, checkIfExists: true)
    referenceGenomeDir = file(referenceGenome.getParent(), checkIfExists: true)
    referenceGenomeIdx = file(params.referenceGenomeIdx, checkIfExists: true)
    referenceGenomeGff = file(params.referenceGenomeGff, checkIfExists: true)
    // databases
    amrfinderDb = file(params.amrfinderDb, checkIfExists: true)
    mlstBlastDb = file(params.mlstBlastDb, checkIfExists: true)
    pubMlstDb = file(params.pubMlstDb, checkIfExists: true)
    chewbbacaDb = file(params.chewbbacaDb, checkIfExists: true)
    coreLociBed = file(params.coreLociBed, checkIfExists: true)
    trainingFile = file(params.trainingFile, checkIfExists: true)
    resfinderDb = file(params.resfinderDb, checkIfExists: true)
    pointfinderDb = file(params.pointfinderDb, checkIfExists: true)
    serotypefinderDb = file(params.serotypefinderDb, checkIfExists: true)
    virulencefinderDb = file(params.virulencefinderDb, checkIfExists: true)

    main:
        ch_versions = Channel.empty()

        CALL_BACTERIAL_BASE( coreLociBed, referenceGenome, referenceGenomeDir, ch_meta.iontorrent, ch_meta.illumina, ch_meta.nanopore )
        
        CALL_BACTERIAL_BASE.out.assembly.set{ch_assembly}
        CALL_BACTERIAL_BASE.out.reads.set{ch_reads}
        CALL_BACTERIAL_BASE.out.bam.set{ch_ref_bam}
        CALL_BACTERIAL_BASE.out.bai.set{ch_ref_bai}
        CALL_BACTERIAL_BASE.out.quast.set{ch_quast}
        CALL_BACTERIAL_BASE.out.qc.set{ch_qc}
        CALL_BACTERIAL_BASE.out.metadata.set{ch_metadata}
        CALL_BACTERIAL_BASE.out.seqrun_meta.set{ch_seqrun_meta}
        CALL_BACTERIAL_BASE.out.input_meta.set{ch_input_meta}
        CALL_BACTERIAL_BASE.out.sourmash.set{ch_sourmash}

        bwa_index(ch_assembly)

        ch_reads
            .join(bwa_index.out.idx)
            .multiMap { id, reads, bai -> 
                reads: tuple(id, reads)
                bai: bai
            }
            .set { bwa_mem_dedup_ch }
        bwa_mem_dedup(bwa_mem_dedup_ch.reads, bwa_mem_dedup_ch.bai)
        samtools_index_assembly(bwa_mem_dedup.out.bam)

        // construct freebayes input channels
        ch_assembly
            .join(bwa_mem_dedup.out.bam)
            .join(samtools_index_assembly.out.bai)
            .multiMap { id, fasta, bam, bai -> 
                assembly: tuple(id, fasta)
                mapping: tuple(bam, bai)
            }
            .set { freebayes_ch }

        // VARIANT CALLING
        freebayes(freebayes_ch.assembly, freebayes_ch.mapping)

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
        serotypefinder(ch_reads, params.useSerotypeDbs, serotypefinderDb)

        // SCREENING
        // antimicrobial detection (amrfinderplus)
        amrfinderplus(ch_assembly, params.species, amrfinderDb)

        // resistance & virulence prediction
        resfinder(ch_reads, params.species, resfinderDb, pointfinderDb)
        virulencefinder(ch_reads, params.useVirulenceDbs, virulencefinderDb)

        ch_reads.map { sampleID, reads -> [ sampleID, [] ] }.set{ ch_empty }

        ch_quast
            .join(ch_qc)
            .join(mlst.out.json)
            .join(chewbbaca_split_results.out.output)
            .join(amrfinderplus.out.output)
            .join(resfinder.out.json)
            .join(resfinder.out.meta)
            .join(serotypefinder.out.json)
            .join(serotypefinder.out.meta)
            .join(virulencefinder.out.json)
            .join(virulencefinder.out.meta)
            .join(ch_empty)
            .join(ch_empty)
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
            create_analysis_result(combinedOutput, referenceGenome, referenceGenomeIdx, referenceGenomeGff)
            ch_versions = ch_versions.mix(kraken.out.versions)
            ch_versions = ch_versions.mix(bracken.out.versions)
        } else {
            combinedOutput.join(ch_empty).set{ combinedOutput }
            create_analysis_result(combinedOutput, referenceGenome, referenceGenomeIdx, referenceGenomeGff)
        }

        create_yaml(create_analysis_result.out.json.join(ch_sourmash), params.speciesDir)

        ch_quast
            .join(ch_qc)
            .join(chewbbaca_split_results.out.output)
            .set{ cdmInput }

        create_cdm_input(cdmInput)

        export_to_cdm(create_cdm_input.out.json.join(ch_seqrun_meta), params.speciesDir)

        copy_to_cron(create_yaml.out.yaml.join(export_to_cdm.out.cdm))

        ch_versions = ch_versions.mix(CALL_BACTERIAL_BASE.out.versions)
        ch_versions = ch_versions.mix(amrfinderplus.out.versions)
        ch_versions = ch_versions.mix(bwa_index.out.versions)
        ch_versions = ch_versions.mix(bwa_mem_dedup.out.versions)
        ch_versions = ch_versions.mix(chewbbaca_allelecall.out.versions)
        ch_versions = ch_versions.mix(create_analysis_result.out.versions)
        ch_versions = ch_versions.mix(freebayes.out.versions)
        ch_versions = ch_versions.mix(mlst.out.versions)
        ch_versions = ch_versions.mix(resfinder.out.versions)
        ch_versions = ch_versions.mix(serotypefinder.out.versions)
        ch_versions = ch_versions.mix(samtools_index_assembly.out.versions)
        ch_versions = ch_versions.mix(virulencefinder.out.versions)

    emit: 
        pipeline_result = create_analysis_result.out.json
        cdm             = export_to_cdm.out.cdm
        cron_yaml       = copy_to_cron.out.yaml
        cron_cdm        = copy_to_cron.out.cdm
        versions        = ch_versions
}
