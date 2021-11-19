#!/usr/bin/env nextflow

// Example of an bacterial analysis pipeline
nextflow.enable.dsl=2

// define workflow parameters
params.outdir = 'test_output'
params.specie = 'saureus'

include { samtools_sort as samtools_sort_one; samtools_sort as samtools_sort_two; samtools_index as samtools_index_one; samtools_index as samtools_index_two } from './nextflow-modules/modules/samtools/main.nf'
include { sambamba_markdup } from './nextflow-modules/modules/sambamba/main.nf'
include { freebayes } from './nextflow-modules/modules/freebayes/main.nf' addParams( args: ['-C', '2', '-F', '0.2', '--pooled-continuous'] )
include { spades } from './nextflow-modules/modules/spades/main.nf' addParams( args: ['--iontorrent', '--only-assembler'] )
include { mask_polymorph_assembly; export_to_cdm; export_to_cgviz } from './nextflow-modules/modules/cmd/main.nf'
include { quast } from './nextflow-modules/modules/quast/main.nf' addParams( args: [] )
include { mlst } from './nextflow-modules/modules/mlst/main.nf' addParams( args: [] )
include { ariba_run } from './nextflow-modules/modules/ariba/main.nf' addParams( outdir: 'ariba_test_outdir', args: ['--force'] )
include { ariba_summary } from './nextflow-modules/modules/ariba/main.nf' addParams( args: ['--col_filter', 'n', '--row_filter', 'n'] )
include { kraken } from './nextflow-modules/modules/kraken/main.nf' addParams( args: ['--gzip-compressed'] )
include { bracken } from './nextflow-modules/modules/bracken/main.nf' addParams( args: ['-r', '150'] )
include { bwa_mem as bwa_mem_ref; bwa_mem as bwa_mem_dedup; bwa_index } from './nextflow-modules/modules/bwa/main.nf' addParams( args: ['-M'] )
include { chewbbaca_allelecall; chewbbaca_split_results; chewbbaca_split_missing_loci } from './nextflow-modules/modules/chewbbaca/main.nf'


process ariba_summary_to_json {
  tag "${summary.simpleName}"
  label "process_low"
  publishDir "${params.outdir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    path summary
    path report
    path reference

  output:
    path("${output}"), emit: output

  script:
    output = "${summary.simpleName}_export.json"
    """
    ariba2json.pl ${reference} ${summary} ${report} > ${output}
    """
}

process post_align_qc {
  tag "${bam.simpleName}"
  label "process_low"
  publishDir "${params.outdir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    path bam
    path index
    path cgmlst

  output:
    path report

  script:
    id = "${bam.simpleName}"
    output = "${id}_bwa.qc"
    """
    postaln_qc.pl ${bam} ${reference} ${id} ${task.cpus} > ${output}
    """
}

workflow bacterial_pipeline {
  reads = Channel .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row -> tuple(row.id, tuple(file(row.read1), file(row.read2))) }

  // load references 
  genomeReference = file(params.genomeReference, checkIfExists: true)
  genomeReferenceDir = file(genomeReference.getParent(), checkIfExists: true)
  aribaReference = file(params.aribaReference, checkIfExists: true)
  aribaReferenceDir = file(aribaReference.getParent(), checkIfExists: true)
  // databases
  krakenDb = file(params.krakenDb, checkIfExists: true)
  cgmlstDb = file(params.cgmlstDb, checkIfExists: true)
  trainingFile = file(params.trainingFile, checkIfExists: true)


  main:
    // assembly and qc processing
    referenceMapping = bwa_mem_ref(reads, genomeReferenceDir)
    sortedReferenceMapping = samtools_sort_one(referenceMapping, [])
    sortedReferenceMappingIdx = samtools_index_one(sortedReferenceMapping.bam)

    sambamba_markdup(sortedReferenceMapping.bam, sortedReferenceMappingIdx)
    postQc = post_align_qc(sortedReferenceMapping, sortedReferenceMappingIdx)

    assembly = spades(reads)

    // mask polymorph regions
    assemblyBwaIdx = bwa_index(assembly)
    assemblyMapping = bwa_mem_dedup(reads, assemblyBwaIdx)
    sortedAssemblyMapping = samtools_sort_two(assemblyMapping, [])
    sortedAssemblyMappingIdx = samtools_index_two(sortedAssemblyMapping.bam)
    maskedRegionsVcf = freebayes(assembly, sortedAssemblyMappingIdx, sortedAssemblyMapping.bam)
    maskedAssembly = mask_polymorph_assembly(assembly, maskedRegionsVcf)

    // typing path
    assemblyQc = quast(assembly, genomeReference)
    mlstResult = mlst(assembly, params.specie)
    chewbbacaResult = chewbbaca_allelecall(maskedAssembly, cgmlstDb, trainingFile)
    chewbbaca_split_results(chewbbacaResult.results)
    chewbbaca_split_missing_loci(chewbbacaResult.missing)

    // end point
    export_to_cdm(chewbbacaResult.missing, assemblyQc, postQc)

    // ariba path
    aribaReport = ariba_run(reads, aribaReferenceDir)
    aribaSummary = ariba_summary(aribaReport)
    aribaJson = ariba_summary_to_json(aribaReport, aribaSummary, aribaReference)

    // kraken path
    krakenReport  = kraken(reads, krakenDb).report
    brackenOutput = bracken(krakenReport, krakenDb).output

    export_to_cgviz(assemblyQc, mlstResult.json, chewbbacaResult.results, krakenReport, aribaJson)
}
