#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { get_meta                              } from '../methods/get_sample_data.nf'
include { get_reads                             } from '../methods/get_sample_data.nf'
include { get_seqrun_meta                       } from '../methods/get_seqrun_meta.nf'
include { assembly_trim_clean                   } from '../nextflow-modules/modules/clean/main.nf'
include { bwa_mem as bwa_mem_ref                } from '../nextflow-modules/modules/bwa/main.nf'
include { fastqc                                } from '../nextflow-modules/modules/fastqc/main.nf'
include { flye                                  } from '../nextflow-modules/modules/flye/main.nf'
include { hostile                               } from '../nextflow-modules/modules/hostile/main.nf'
include { medaka                                } from '../nextflow-modules/modules/medaka/main.nf'
include { nanoplot                              } from '../nextflow-modules/modules/nanoplot/main.nf'
include { minimap2_to_ref                       } from '../nextflow-modules/modules/minimap2/main.nf'       
include { post_align_qc                         } from '../nextflow-modules/modules/prp/main.nf'
include { quast                                 } from '../nextflow-modules/modules/quast/main.nf'
include { samtools_index as samtools_index_ref  } from '../nextflow-modules/modules/samtools/main.nf'
include { samtools_coverage                     } from '../nextflow-modules/modules/samtools/main.nf'
include { samtools_sort                         } from '../nextflow-modules/modules/samtools/main.nf'
include { save_analysis_metadata                } from '../nextflow-modules/modules/meta/main.nf'
include { seqtk_sample                          } from '../nextflow-modules/modules/seqtk/main.nf'
include { ska_build                             } from '../nextflow-modules/modules/ska/main.nf'
include { skesa                                 } from '../nextflow-modules/modules/skesa/main.nf'
include { sourmash                              } from '../nextflow-modules/modules/sourmash/main.nf'
include { spades_illumina                       } from '../nextflow-modules/modules/spades/main.nf'
include { spades_iontorrent                     } from '../nextflow-modules/modules/spades/main.nf'

workflow CALL_BACTERIAL_BASE {
    take:
        coreLociBed
        referenceGenome
        referenceGenomeDir
        referenceGenomeMmi
        inputSamples
        targetSampleSize
    
    main:
        ch_versions = Channel.empty()

        // Create channel for sample metadata
        Channel.fromPath(inputSamples)
            .splitCsv(header:true)
            .tap{ ch_raw_input }
            .map{ row -> get_meta(row) }
            .set{ ch_meta }

        // Create channel for reads
        ch_raw_input
            .map{ row -> get_reads(row) }
            .set{ ch_raw_reads }

        if ( params.useHostile ) {
            // remove human reads
            hostile( ch_raw_reads ).reads.set{ ch_depleted_reads }
            ch_versions = ch_versions.mix(hostile.out.versions)
        } else {
            ch_raw_reads.set{ ch_depleted_reads }
        }

        if ( params.targetSampleSize ) {
            // downsample reads
            seqtk_sample( ch_depleted_reads, targetSampleSize ).reads.set{ ch_depleted_sampled_reads }
            ch_versions = ch_versions.mix(seqtk_sample.out.versions)
        } else {
            ch_depleted_reads.set{ ch_depleted_sampled_reads }
        }

        // reads trim and clean and recreate reads channel if the reads were filtered or downsampled
        assembly_trim_clean(ch_depleted_sampled_reads.join(ch_meta)).set { ch_clean_reads_w_meta }
        Channel.empty()
            .mix( ch_depleted_sampled_reads, ch_clean_reads_w_meta )  // if samples are filtered or downsampled
            .tap{ ch_reads }                                          // create reads channel
            .join( ch_meta )                                          // add meta info
            .set{ ch_reads_w_meta }                                   // write as temp channel

        if ( params.cronCopy || params.devMode) {
            Channel.fromPath(params.csv).splitCsv(header:true)
                .map{ row -> get_seqrun_meta(row) }
                .set{ ch_seqrun_meta }
        } else {
            ch_reads.map{ sampleID, reads -> [ sampleID, [], [], [] ] }.set{ ch_seqrun_meta }
        }
        // analysis metadata
        save_analysis_metadata(ch_reads_w_meta.join(ch_seqrun_meta))

        // assembly
        skesa(ch_reads_w_meta)
        spades_illumina(ch_reads_w_meta)
        spades_iontorrent(ch_reads_w_meta)
        flye(ch_reads_w_meta)
        medaka(ch_reads_w_meta, flye.out.fasta)

        Channel.empty()
            .mix(
                skesa.out.fasta, spades_illumina.out.fasta, 
                spades_iontorrent.out.fasta, medaka.out.fasta
            ).set{ ch_assembly }

        // evaluate assembly quality 
        quast(ch_assembly, referenceGenome)

        // qc processing - short read
        fastqc(ch_reads)
        bwa_mem_ref(ch_reads, referenceGenomeDir)
        samtools_index_ref(bwa_mem_ref.out.bam)

        post_align_qc(bwa_mem_ref.out.bam, referenceGenome, coreLociBed)

        // qc processing - long read
        nanoplot(ch_reads_w_meta)
        minimap2_to_ref(ch_reads_w_meta, referenceGenomeMmi)
        samtools_sort(minimap2_to_ref.out.sam)
        //index the sorted bam file
        //samtools_index_ref(samtools_sort.out.bam) - but then it gets in the same folder as SR bam, which needs to be avoided by blocking bwa_mem_ref for analysis of LR
        samtools_coverage(samtools_sort.out.bam)


        sourmash(ch_assembly)

        ska_build(ch_reads)

        ch_versions = ch_versions.mix(bwa_mem_ref.out.versions)
        ch_versions = ch_versions.mix(fastqc.out.versions)
        ch_versions = ch_versions.mix(flye.out.versions)
        ch_versions = ch_versions.mix(medaka.out.versions)
        ch_versions = ch_versions.mix(minimap2_to_ref.out.versions)
        ch_versions = ch_versions.mix(nanoplot.out.versions)
        ch_versions = ch_versions.mix(quast.out.versions)
        ch_versions = ch_versions.mix(samtools_index_ref.out.versions)
        ch_versions = ch_versions.mix(skesa.out.versions)
        ch_versions = ch_versions.mix(ska_build.out.versions)
        ch_versions = ch_versions.mix(sourmash.out.versions)
        ch_versions = ch_versions.mix(spades_illumina.out.versions)
        ch_versions = ch_versions.mix(spades_iontorrent.out.versions)

    emit:
        assembly        = ch_assembly                       // channel: [ val(meta), path(fasta)]
        bam             = bwa_mem_ref.out.bam               // channel: [ val(meta), path(bam)]
        bai             = samtools_index_ref.out.bai        // channel: [ val(meta), path(bai)]
        fastqc          = fastqc.out.output                 // channel: [ val(meta), path(txt)]
        seqplat_meta    = ch_meta                           // channel: [ val(meta), val(str)]
        metadata        = save_analysis_metadata.out.meta   // channel: [ val(meta), path(json)]
        qc              = post_align_qc.out.qc              // channel: [ val(meta), path(fasta)]
        qc_nano_raw     = nanoplot.out.html                 // channel: [ val(meta), path(html)]
        qc_nano_cov     = samtools_coverage.out.txt         // channel: [ val(meta), path(txt)]
        quast           = quast.out.qc                      // channel: [ val(meta), path(qc)]
        reads           = ch_reads                          // channel: [ val(meta), path(json)]
        reads_w_meta    = ch_reads_w_meta                   // channel: [ val(meta), path(meta)]
        ska_build       = ska_build.out.skf                 // channel: [ val(meta), path(skf)]
        seqrun_meta     = ch_seqrun_meta                    // channel: [ val(meta), val(json), val(json)]
        sourmash        = sourmash.out.signature            // channel: [ val(meta), path(signature)]
        versions        = ch_versions                       // channel: [ versions.yml ]
}
