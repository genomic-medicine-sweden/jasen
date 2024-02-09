include { chewbbaca_allelecall          } from '../nextflow-modules/modules/chewbbaca/main.nf'
include { chewbbaca_create_batch_list   } from '../nextflow-modules/modules/chewbbaca/main.nf'
include { chewbbaca_split_results       } from '../nextflow-modules/modules/chewbbaca/main.nf'
include { mlst                          } from '../nextflow-modules/modules/mlst/main.nf'
include { mask_polymorph_assembly       } from '../nextflow-modules/modules/mask/main.nf'

workflow CALL_TYPING {
    take:
        chewbbacaDb
        mlstDb
        trainingFile
        ch_assembly         // channel: [ val(meta), val(fasta) ]
        ch_freebayes_vcf    // channel: [ val(meta), val(vcf) ]

    main:
        ch_versions = Channel.empty()

        mask_polymorph_assembly(ch_assembly.join(ch_freebayes_vcf))
        // sequence typing
        mlst(ch_assembly, params.species, mlstDb)
        // split assemblies and id into two seperate channels to enable re-pairing
        // of results and id at a later stage. This to allow batch cgmlst/wgmlst analysis 
        mask_polymorph_assembly.out.fasta
            .multiMap { sampleName, filePath -> 
                sampleName: sampleName
                filePath: filePath
            }
            .set{ maskedAssemblyMap }

        chewbbaca_create_batch_list(maskedAssemblyMap.filePath.collect())
        chewbbaca_allelecall(maskedAssemblyMap.sampleName.collect(), chewbbaca_create_batch_list.out.list, chewbbacaDb, trainingFile)
        chewbbaca_split_results(chewbbaca_allelecall.out.sampleName, chewbbaca_allelecall.out.calls)

        ch_versions = ch_versions.mix(mlst.out.versions)
        ch_versions = ch_versions.mix(chewbbaca_allelecall.out.versions)

    emit:
        chewbbaca = chewbbaca_split_results.out.output  // channel: [ val(meta), path(tsv)]
        mlst      = mlst.out.json                       // channel: [ val(meta), path(json)]
        versions  = ch_versions                         // channel: [ versions.yml ]
}