#!/usr/bin/env nextflow

// Example of an bacterial analysis pipeline
nextflow.enable.dsl=2

include { samtools_sort as samtools_sort_one; samtools_sort as samtools_sort_two; samtools_index as samtools_index_one; samtools_index as samtools_index_two } from './nextflow-modules/modules/samtools/main.nf'
include { sambamba_markdup } from './nextflow-modules/modules/sambamba/main.nf'
include { freebayes } from './nextflow-modules/modules/freebayes/main.nf' addParams( args: ['-C', '2', '-F', '0.2', '--pooled-continuous'] )
include { spades } from './nextflow-modules/modules/spades/main.nf' addParams( args: ['--only-assembler'] )
include { skesa } from './nextflow-modules/modules/skesa/main.nf'
include { save_analysis_metadata; mask_polymorph_assembly; export_to_cdm } from './nextflow-modules/modules/cmd/main.nf'
include { quast } from './nextflow-modules/modules/quast/main.nf' addParams( args: [] )
include { mlst } from './nextflow-modules/modules/mlst/main.nf' addParams( args: [] )
include { ariba_prepareref } from './nextflow-modules/modules/ariba/main.nf'
include { ariba_run } from './nextflow-modules/modules/ariba/main.nf' addParams( args: ['--force'] )
include { ariba_summary } from './nextflow-modules/modules/ariba/main.nf' addParams( args: ['--col_filter', 'n', '--row_filter', 'n'] )
include { ariba_summary_to_json } from './nextflow-modules/modules/ariba/main.nf'
include { kraken } from './nextflow-modules/modules/kraken/main.nf' addParams( args: ['--gzip-compressed'] )
include { bracken } from './nextflow-modules/modules/bracken/main.nf' addParams( args: ['-r', '150'] )
include { bwa_mem as bwa_mem_ref; bwa_mem as bwa_mem_dedup; bwa_index } from './nextflow-modules/modules/bwa/main.nf' addParams( args: ['-M'] )
include { post_align_qc } from './nextflow-modules/modules/qc/main.nf'
include { chewbbaca_allelecall; chewbbaca_split_results; chewbbaca_split_missing_loci; chewbbaca_create_batch_list} from './nextflow-modules/modules/chewbbaca/main.nf' addParams( args: [] )
include { resfinder } from './nextflow-modules/modules/resfinder/main.nf' addParams( args: [] )
include { virulencefinder } from './nextflow-modules/modules/virulencefinder/main.nf' addParams( args: [] )
include { create_analysis_result } from './nextflow-modules/modules/prp/main.nf' addParams( args: [] )

workflow bacterial_default {
  if (params.strands==2){
    reads = Channel.fromPath(params.csv)
    .splitCsv(header:true).map{ row -> tuple(row.id, tuple(file(row.read1),file(row.read2))) }}
  else if (params.strands==1){
    reads = Channel.fromPath(params.csv)
    .splitCsv(header:true).map{ row -> tuple(row.id, tuple(file(row.read1))) }}

  // load references 
  genomeReference = file(params.genomeReference, checkIfExists: true)
  genomeReferenceDir = file(genomeReference.getParent(), checkIfExists: true)
  aribaReference = file(params.aribaReference, checkIfExists: true)
  aribaReferenceDir = file(aribaReference.getParent(), checkIfExists: true)
  // databases
  mlstDb = file(params.mlstBlastDb, checkIfExists: true)
  cgmlstDb = file(params.cgmlstDb, checkIfExists: true)
  cgmlstSchema = file(params.cgmlstSchema, checkIfExists: true)
  trainingFile = file(params.trainingFile, checkIfExists: true)
  resfinderDb = file(params.resfinderDb, checkIfExists: true)
  pointfinderDb = file(params.pointfinderDb, checkIfExists: true)
  virulencefinderDb = file(params.virulencefinderDb, checkIfExists: true)


  main:
    runInfo = save_analysis_metadata()
    // assembly and qc processing
    referenceMapping = bwa_mem_ref(reads, genomeReferenceDir)
    sortedReferenceMapping = samtools_sort_one(referenceMapping, [])
    sortedReferenceMappingIdx = samtools_index_one(sortedReferenceMapping.bam)

    //sambamba_markdup(sortedReferenceMapping.bam, sortedReferenceMappingIdx)
    sortedReferenceMapping.bam
      .join(sortedReferenceMappingIdx)
      .multiMap { id, bam, bai -> 
          bam: tuple(id, bam)
          bai: bai
       }
      .set{ post_align_qc_ch }
    postQc = post_align_qc(post_align_qc_ch.bam, post_align_qc_ch.bai, cgmlstSchema)
    
    if ( params.useSkesa ) {
      assembly = skesa(reads)
    } else {
      assembly = spades(reads)
    }

    // mask polymorph regions
    assemblyBwaIdx = bwa_index(assembly)
    reads
      .join(assemblyBwaIdx)
      .multiMap { id, reads, bai -> 
          reads: tuple(id, reads)
          bai: bai
       }
      .set { bwa_mem_dedup_ch }
    assemblyMapping = bwa_mem_dedup(bwa_mem_dedup_ch.reads, bwa_mem_dedup_ch.bai)
    sortedAssemblyMapping = samtools_sort_two(assemblyMapping, [])
    sortedAssemblyMappingIdx = samtools_index_two(sortedAssemblyMapping.bam)
    // construct freebayes input channels
    assembly
      .join(sortedAssemblyMapping.bam)
      .join(sortedAssemblyMappingIdx)
      .multiMap { id, fasta, bam, bai -> 
          assembly: tuple(id, fasta)
          mapping: tuple(bam, bai)
       }
      .set { freebayes_ch }
    maskedRegionsVcf = freebayes(freebayes_ch.assembly, freebayes_ch.mapping)
    //maskedRegionsVcf = freebayes(assembly.join(sortedAssemblyMappingIdx).join(sortedAssemblyMapping.bam))
    maskedAssembly = mask_polymorph_assembly(assembly.join(maskedRegionsVcf))
    // maskedAssemblies = maskedAssembly.map({ sampleName, filePath -> [ filePath ] }).collect()
    assemblyQc = quast(assembly, genomeReference)
    mlstResult = mlst(assembly, params.specie, mlstDb)
    // split assemblies and id into two seperate channels to enable re-pairing
    // of results and id at a later stage. This to allow batch cgmlst analysis 
    maskedAssembly
      .multiMap { sampleName, filePath -> 
          sampleName: sampleName
          filePath: filePath
      }
      .set{ maskedAssemblyMap }

    batchList = chewbbaca_create_batch_list(maskedAssemblyMap.filePath.collect())
    chewbbaca_allelecall(maskedAssemblyMap.sampleName.collect(), batchList, cgmlstDb, trainingFile)
    chewbbacaResult = chewbbaca_split_results(chewbbaca_allelecall.out.sampleName, chewbbaca_allelecall.out.calls)
    //chewbbaca_split_missing_loci(chewbbacaResult.missing)

    // end point
    export_to_cdm(chewbbacaResult.join(assemblyQc).join(postQc))

    // ariba path
    if (params.strands == 2) {
      // debug aribaReferenceDir = ariba_prepareref(aribaReference, Channel.empty())
      aribaReport = ariba_run(reads, aribaReferenceDir)
      aribaSummary = ariba_summary(aribaReport)
      aribaJson = ariba_summary_to_json(aribaReport.join(aribaSummary), aribaReference)
    }

    // perform resistance prediction
    resfinderOutput = resfinder(reads, params.specie, resfinderDb, pointfinderDb)
    virulencefinderOutput = virulencefinder(reads, params.useVirulenceDbs, virulencefinderDb)

    // combine results for export
    combinedOutput = assemblyQc
        .join(mlstResult.json)
        .join(chewbbacaResult)
        .join(resfinderOutput.json)
        .join(resfinderOutput.meta)
        .join(virulencefinderOutput.json)
        .join(virulencefinderOutput.meta)

    // Using kraken for species identificaiton
    if( params.useKraken ) {
      krakenDb = file(params.krakenDb, checkIfExists: true)
      krakenReport = kraken(reads, krakenDb).report
      brackenOutput = bracken(krakenReport, krakenDb).output
      combinedOutput = combinedOutput.join(brackenOutput)
      create_analysis_result(
        runInfo, 
        combinedOutput
      )
	  } else {
      emptyBrackenOutput = maskedAssemblyMap.sampleName.map{sampleName -> [ sampleName, [] ] }
      combinedOutput = combinedOutput.join(emptyBrackenOutput)
      create_analysis_result(
        runInfo, 
        combinedOutput
      )
	  }
    

  emit: 
    pipeline_result = create_analysis_result.output
    cdm_import = export_to_cdm.output
}
