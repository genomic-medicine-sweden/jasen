#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { get_meta              } from '../methods/get_meta'
include { CALL_ASSEMBLY         } from '../subworkflows/assembly.nf'
include { CALL_POSTPROCESSING   } from '../subworkflows/postprocessing.nf'
include { CALL_PREPROCESSING    } from '../subworkflows/preprocessing.nf'
include { CALL_QUALITY_CONTROL  } from '../subworkflows/quality_control.nf'
include { CALL_RELATEDNESS      } from '../subworkflows/relatedness.nf'
include { CALL_SCREENING        } from '../subworkflows/screening.nf'
include { CALL_TYPING           } from '../subworkflows/typing.nf'
include { CALL_VARIANT_CALLING  } from '../subworkflows/variant_calling.nf'

workflow CALL_BACTERIAL_DEFAULT {
    Channel.fromPath(params.csv).splitCsv(header:true)
        .map{ row -> get_meta(row) }
        .branch {
        iontorrent: it[2] == "iontorrent"
        illumina: it[2] == "illumina"
        }
        .set{ ch_meta }

    // load references 
    genomeReference = file(params.genomeReference, checkIfExists: true)
    genomeReferenceDir = file(genomeReference.getParent(), checkIfExists: true)
    // databases
    amrfinderDb = file(params.amrfinderDb, checkIfExists: true)
    mlstDb = file(params.mlstBlastDb, checkIfExists: true)
    chewbbacaDb = file(params.chewbbacaDb, checkIfExists: true)
    coreLociBed = file(params.coreLociBed, checkIfExists: true)
    trainingFile = file(params.trainingFile, checkIfExists: true)
    resfinderDb = file(params.resfinderDb, checkIfExists: true)
    pointfinderDb = file(params.pointfinderDb, checkIfExists: true)
    virulencefinderDb = file(params.virulencefinderDb, checkIfExists: true)

    main:
        CALL_PREPROCESSING( ch_meta.iontorrent, ch_meta.illumina )
        CALL_ASSEMBLY( genomeReference, CALL_PREPROCESSING.out.input_meta )
        CALL_VARIANT_CALLING( CALL_ASSEMBLY.out.assembly, CALL_PREPROCESSING.out.input_meta, CALL_PREPROCESSING.out.reads )
        CALL_QUALITY_CONTROL( coreLociBed, genomeReferenceDir, CALL_PREPROCESSING.out.reads )
        CALL_TYPING( chewbbacaDb, mlstDb, trainingFile, CALL_ASSEMBLY.out.assembly, CALL_VARIANT_CALLING.out.vcf )
        CALL_SCREENING( amrfinderDb, pointfinderDb, resfinderDb, virulencefinderDb, CALL_ASSEMBLY.out.assembly, CALL_PREPROCESSING.out.reads )
        CALL_RELATEDNESS( CALL_ASSEMBLY.out.assembly )

        CALL_ASSEMBLY.out.qc
            .join(CALL_QUALITY_CONTROL.out.qc)
            .join(CALL_TYPING.out.mlst)
            .join(CALL_TYPING.out.chewbbaca)
            .join(CALL_SCREENING.out.amrfinderplus)
            .join(CALL_SCREENING.out.resfinderJson)
            .join(CALL_SCREENING.out.resfinderMeta)
            .join(CALL_SCREENING.out.virulencefinderJson)
            .join(CALL_SCREENING.out.virulencefinderMeta)
            .join(CALL_PREPROCESSING.out.metadata)
            .set{ combinedOutput }

        CALL_POSTPROCESSING( CALL_TYPING.out.chewbbaca, combinedOutput, CALL_QUALITY_CONTROL.out.qc, CALL_ASSEMBLY.out.qc, CALL_PREPROCESSING.out.reads )

        emit: 
            pipeline_result = CALL_POSTPROCESSING.out.pipeline_result
            cdm             = CALL_POSTPROCESSING.out.cdm
}