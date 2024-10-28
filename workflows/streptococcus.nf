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
include { create_analysis_result                    } from '../nextflow-modules/modules/prp/main.nf'
include { create_cdm_input                          } from '../nextflow-modules/modules/prp/main.nf'
include { create_yaml                               } from '../nextflow-modules/modules/yaml/main.nf'
include { emmtyper                                  } from '../nextflow-modules/modules/emmtyper/main.nf'
include { export_to_cdm                             } from '../nextflow-modules/modules/cmd/main.nf'
include { freebayes                                 } from '../nextflow-modules/modules/freebayes/main.nf'
include { kraken                                    } from '../nextflow-modules/modules/kraken/main.nf'
include { mask_polymorph_assembly                   } from '../nextflow-modules/modules/mask/main.nf'
include { mlst                                      } from '../nextflow-modules/modules/mlst/main.nf'
include { resfinder                                 } from '../nextflow-modules/modules/resfinder/main.nf'
include { samtools_index as samtools_index_assembly } from '../nextflow-modules/modules/samtools/main.nf'
include { serotypefinder                            } from '../nextflow-modules/modules/serotypefinder/main.nf'
include { shigapass                                 } from '../nextflow-modules/modules/shigapass/main.nf'
include { virulencefinder                           } from '../nextflow-modules/modules/virulencefinder/main.nf'
include { CALL_BACTERIAL_BASE                       } from '../workflows/bacterial_base.nf'

workflow CALL_STREPTOCOCCUS {
    Channel.fromPath(params.csv)
        .splitCsv(header:true)
        .map{ row -> get_meta(row) }
        .branch {
            iontorrent: it[2] == "iontorrent"
            illumina: it[2] == "illumina"
            nanopore: it[2] == "nanopore"
        }
        .set{ ch_meta }

    // load references
    referenceGenome = params.referenceGenome ? file(params.referenceGenome, checkIfExists: true) : Channel.of([])
    referenceGenomeDir = params.referenceGenome ? file(referenceGenome.getParent(), checkIfExists: true) : Channel.of([])
    referenceGenomeGff = params.referenceGenomeGff ? file(params.referenceGenomeGff, checkIfExists: true) : Channel.of([])
    referenceGenomeIdx = params.referenceGenomeIdx ? file(params.referenceGenomeIdx, checkIfExists: true) : Channel.of([])
    // databases
    amrfinderDb = file(params.amrfinderDb, checkIfExists: true)
    chewbbacaDb = file(params.chewbbacaDb, checkIfExists: true)
    coreLociBed = params.coreLociBed ? file(params.coreLociBed, checkIfExists: true) : Channel.of([])
    krakenDb = params.krakenDb ? file(params.krakenDb, checkIfExists: true) : Channel.of([])
    mlstBlastDb = params.mlstBlastDb ? file(params.mlstBlastDb, checkIfExists: true) : Channel.of([])
    pointfinderDb = file(params.pointfinderDb, checkIfExists: true)
    pubMlstDb = params.pubMlstDb ? file(params.pubMlstDb, checkIfExists: true) : Channel.of([])
    resfinderDb = file(params.resfinderDb, checkIfExists: true)
    serotypefinderDb = file(params.serotypefinderDb, checkIfExists: true)
    shigapassDb = params.shigapassDb ? file(params.shigapassDb, checkIfExists: true) : Channel.of([])
    trainingFile = params.trainingFile ? file(params.trainingFile, checkIfExists: true) : Channel.of([])
    virulencefinderDb = file(params.virulencefinderDb, checkIfExists: true)
    // schemas and values
    mlstScheme = params.mlstScheme ? params.mlstScheme : Channel.of([])
    species = params.species ? params.species : Channel.of([])
    speciesDir = params.speciesDir ? params.speciesDir : Channel.of([])

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
        CALL_BACTERIAL_BASE.out.ska_build.set{ch_ska}

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
        mlst(ch_assembly, mlstScheme, pubMlstDb, mlstBlastDb)

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
        serotypefinder(ch_reads, params.useSerotypeDbs, serotypefinderDb)
        shigapass(ch_assembly, shigapassDb)

        // SCREENING
        // antimicrobial detection (amrfinderplus)
        amrfinderplus(ch_assembly, species, amrfinderDb)

        // resistance & virulence prediction
        resfinder(ch_reads, species, resfinderDb, pointfinderDb)
        virulencefinder(ch_reads, params.useVirulenceDbs, virulencefinderDb)

        ch_reads.map { sampleID, reads -> [ sampleID, [] ] }.set{ ch_empty }

        ch_amrfinderplus = amrfinderplus.out.output.concat(ch_empty).first()
        ch_chewbbaca = chewbbaca_split_results.out.output.concat(ch_empty).first()
        ch_emmtyper = emmtyper.out.tsv.concat(ch_empty).first()
        ch_mlst = mlst.out.json.concat(ch_empty).first()
        ch_quast = ch_quast.concat(ch_empty).first()
        ch_qc = ch_qc.concat(ch_empty).first()
        ch_ref_bam = ch_ref_bam.concat(ch_empty).first()
        ch_ref_bai = ch_ref_bai.concat(ch_empty).first()
        ch_resfinder = resfinder.out.json.concat(ch_empty).first()
        ch_resfinder_meta = resfinder.out.meta.concat(ch_empty).first()
        ch_serotypefinder = serotypefinder.out.json.concat(ch_empty).first()
        ch_serotypefinder_meta = serotypefinder.out.meta.concat(ch_empty).first()
        ch_shigapass = shigapass.out.csv.concat(ch_empty).first()
        ch_virulencefinder = virulencefinder.out.json.concat(ch_empty).first()
        ch_virulencefinder_meta = virulencefinder.out.meta.concat(ch_empty).first()

        ch_quast
            .join(ch_qc)
            .join(ch_mlst)
            .join(ch_chewbbaca)
            .join(ch_amrfinderplus)
            .join(ch_resfinder)
            .join(ch_resfinder_meta)
            .join(ch_serotypefinder)
            .join(ch_serotypefinder_meta)
            .join(ch_virulencefinder)
            .join(ch_virulencefinder_meta)
            .join(ch_shigapass)
            .join(ch_emmtyper)
            .join(ch_ref_bam)
            .join(ch_ref_bai)
            .join(ch_metadata)
            .join(ch_empty)
            .join(ch_empty)
            .join(ch_empty)
            .set{ combinedOutput }

        kraken(ch_reads, krakenDb)
        bracken(kraken.out.report, krakenDb)
        brackenOutput = bracken.out.output ? bracken.out.output : ch_empty

        combinedOutput.join(bracken.out.output).set{ combinedOutput }
        create_analysis_result(combinedOutput, referenceGenome, referenceGenomeIdx, referenceGenomeGff)

        create_yaml(create_analysis_result.out.json.join(ch_sourmash).join(ch_ska), speciesDir)

        ch_quast
            .join(ch_qc)
            .join(chewbbaca_split_results.out.output)
            .set{ cdmInput }

        create_cdm_input(cdmInput)

        export_to_cdm(create_cdm_input.out.json.join(ch_seqrun_meta), speciesDir)

        ch_versions = ch_versions.mix(CALL_BACTERIAL_BASE.out.versions)
        ch_versions = ch_versions.mix(amrfinderplus.out.versions)
        ch_versions = ch_versions.mix(bracken.out.versions)
        ch_versions = ch_versions.mix(bwa_index.out.versions)
        ch_versions = ch_versions.mix(bwa_mem_dedup.out.versions)
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
