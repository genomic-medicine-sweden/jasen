#!/usr/bin/env nextflow

nextflow.enable.dsl=2
idList = []     // Global variable to save IDs from samplesheet

include { abritamr                                  } from './nextflow-modules/modules/abritamr/main.nf'
include { amrfinderplus                             } from './nextflow-modules/modules/amrfinderplus/main.nf'
include { assembly_trim_clean                       } from './nextflow-modules/modules/clean/main.nf'
include { bracken                                   } from './nextflow-modules/modules/bracken/main.nf'
include { bwa_mem as bwa_mem_ref                    } from './nextflow-modules/modules/bwa/main.nf'
include { bwa_mem as bwa_mem_dedup                  } from './nextflow-modules/modules/bwa/main.nf'
include { bwa_index                                 } from './nextflow-modules/modules/bwa/main.nf'
include { chewbbaca_allelecall                      } from './nextflow-modules/modules/chewbbaca/main.nf'
include { chewbbaca_split_results                   } from './nextflow-modules/modules/chewbbaca/main.nf'
include { chewbbaca_create_batch_list               } from './nextflow-modules/modules/chewbbaca/main.nf'
include { create_analysis_result                    } from './nextflow-modules/modules/prp/main.nf'
include { export_to_cdm                             } from './nextflow-modules/modules/cmd/main.nf'
include { freebayes                                 } from './nextflow-modules/modules/freebayes/main.nf'
include { kraken                                    } from './nextflow-modules/modules/kraken/main.nf'
include { mask_polymorph_assembly                   } from './nextflow-modules/modules/mask/main.nf'
include { mlst                                      } from './nextflow-modules/modules/mlst/main.nf'
include { post_align_qc                             } from './nextflow-modules/modules/qc/main.nf'
include { quast                                     } from './nextflow-modules/modules/quast/main.nf'
include { resfinder                                 } from './nextflow-modules/modules/resfinder/main.nf'
include { samtools_index as samtools_index_ref      } from './nextflow-modules/modules/samtools/main.nf'
include { samtools_index as samtools_index_assembly } from './nextflow-modules/modules/samtools/main.nf'
include { save_analysis_metadata                    } from './nextflow-modules/modules/meta/main.nf'
include { sourmash                                  } from './nextflow-modules/modules/sourmash/main.nf'
include { skesa                                     } from './nextflow-modules/modules/skesa/main.nf'
include { spades_illumina                           } from './nextflow-modules/modules/spades/main.nf'
include { spades_iontorrent                         } from './nextflow-modules/modules/spades/main.nf'
include { virulencefinder                           } from './nextflow-modules/modules/virulencefinder/main.nf'



// Function for platform and paired-end or single-end
def get_meta(LinkedHashMap row) {
  platforms = ["illumina", "nanopore", "pacbio", "iontorrent"]
  idErrorMessage = "ERROR: Please check input samplesheet -> ID '${row.id}' is not a unique name in the samplesheet."
  platformErrorMessage = "ERROR: Please check input samplesheet -> Platform is not one of the following:!\n-${platforms.join('\n-')}"
  
  idList.add(row.id)
  identicalIds=idList.clone().unique().size() != idList.size()
  correctPlatform = (row.platform in platforms)
  
  if (correctPlatform && !identicalIds) {
    if (row.read2) {
      meta = tuple(row.id, tuple(file(row.read1), file(row.read2)), row.platform)
    } else {
      meta = tuple(row.id, tuple(file(row.read1)), row.platform)
    }
  } else if (!correctPlatform && identicalIds) {
    exit 1, platformErrorMessage + "\n\n" + idErrorMessage
  } else if (!correctPlatform && !identicalIds) {
    exit 1, platformErrorMessage
  } else if (correctPlatform && identicalIds) {
    exit 1, idErrorMessage
  }
  return meta
}

workflow bacterial_default {
  Channel.fromPath(params.csv).splitCsv(header:true)
    .map{ row -> get_meta(row) }
    .branch {
      iontorrent: it[2] == "iontorrent"
      illumina: it[2] == "illumina"
    }
    .set{ meta }

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
    // reads trim and clean
    assembly_trim_clean(meta.iontorrent).set{ clean_meta }
    meta.illumina.mix(clean_meta).set{ input_meta }
    input_meta.map { sampleName, reads, platform -> [ sampleName, reads ] }.set{ reads }

    // analysis metadata
    save_analysis_metadata(input_meta)

    // assembly and qc processing
    bwa_mem_ref(reads, genomeReferenceDir)
    samtools_index_ref(bwa_mem_ref.out.bam)

    bwa_mem_ref.out.bam
      .join(samtools_index_ref.out.bai)
      .multiMap { id, bam, bai -> 
        bam: tuple(id, bam)
        bai: bai
      }
      .set{ post_align_qc_ch }
    post_align_qc(post_align_qc_ch.bam, post_align_qc_ch.bai, coreLociBed)
    
    // assembly
    skesa(input_meta)
    spades_illumina(input_meta)
    spades_iontorrent(input_meta)

    Channel.empty().mix(skesa.out.fasta, spades_illumina.out.fasta, spades_iontorrent.out.fasta).set{ assembly }

    // mask polymorph regions
    bwa_index(assembly)
    reads
      .join(bwa_index.out.idx)
      .multiMap { id, reads, bai -> 
        reads: tuple(id, reads)
        bai: bai
      }
      .set { bwa_mem_dedup_ch }
    bwa_mem_dedup(bwa_mem_dedup_ch.reads, bwa_mem_dedup_ch.bai)
    samtools_index_assembly(bwa_mem_dedup.out.bam)

    // construct freebayes input channels
    assembly
      .join(bwa_mem_dedup.out.bam)
      .join(samtools_index_assembly.out.bai)
      .multiMap { id, fasta, bam, bai -> 
        assembly: tuple(id, fasta)
        mapping: tuple(bam, bai)
      }
      .set { freebayes_ch }

    freebayes(freebayes_ch.assembly, freebayes_ch.mapping)
    mask_polymorph_assembly(assembly.join(freebayes.out.vcf))
    quast(assembly, genomeReference)
    mlst(assembly, params.species, mlstDb)
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

    // end point
    export_to_cdm(chewbbaca_split_results.out.output.join(quast.out.qc).join(post_align_qc.out.qc))

    // sourmash
    sourmash(assembly)

    // antimicrobial detection (amrfinderplus & abritamr)
    amrfinderplus(assembly, amrfinderDb)
    //abritamr(amrfinderplus.out.output)

    // perform resistance prediction
    resfinder(reads, params.species, resfinderDb, pointfinderDb)
    virulencefinder(reads, params.useVirulenceDbs, virulencefinderDb)

    // combine results for export
    quast.out.qc
      .join(post_align_qc.out.qc)
      .join(mlst.out.json)
      .join(chewbbaca_split_results.out.output)
      .join(amrfinderplus.out.output)
      .join(resfinder.out.json)
      .join(resfinder.out.meta)
      .join(virulencefinder.out.json)
      .join(virulencefinder.out.meta)
      .join(save_analysis_metadata.out.meta)
      .set{ combinedOutput }

    // Using kraken for species identificaiton
    if ( params.useKraken ) {
      krakenDb = file(params.krakenDb, checkIfExists: true)
      kraken(reads, krakenDb)
      bracken(kraken.out.report, krakenDb).output
      combinedOutput = combinedOutput.join(bracken.out.output)
      create_analysis_result(combinedOutput)
	  } else {
      emptyBrackenOutput = reads.map { sampleName, reads -> [ sampleName, [] ] }
      combinedOutput = combinedOutput.join(emptyBrackenOutput)
      create_analysis_result(combinedOutput)
	  }
    
  emit: 
    pipeline_result = create_analysis_result.output
    cdm_import = export_to_cdm.output
}
