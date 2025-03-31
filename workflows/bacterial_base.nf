#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { get_meta                                   } from '../methods/get_sample_data.nf'
include { get_reads                                  } from '../methods/get_sample_data.nf'
include { get_seqrun_meta                            } from '../methods/get_seqrun_meta.nf'
include { assembly_trim_clean                        } from '../modules/local/clean/main.nf'
include { bracken                                    } from '../modules/nf-core/bracken/main.nf'
include { bwa_mem as bwa_mem_ref                     } from '../modules/nf-core/bwa/main.nf'
include { fastqc                                     } from '../modules/nf-core/fastqc/main.nf'
include { flye                                       } from '../modules/nf-core/flye/main.nf'
include { hostile                                    } from '../modules/nf-core/hostile/main.nf'
include { kraken                                     } from '../modules/nf-core/kraken/main.nf'
include { medaka                                     } from '../modules/nf-core/medaka/main.nf'
include { nanoplot                                   } from '../modules/nf-core/nanoplot/main.nf'
include { minimap2_align as minimap2_align_ref       } from '../modules/nf-core/minimap2/main.nf'       
include { post_align_qc                              } from '../modules/local/prp/main.nf'
include { quast                                      } from '../modules/nf-core/quast/main.nf'
include { samtools_index as samtools_index_ref       } from '../modules/nf-core/samtools/main.nf'
include { samtools_coverage as samtools_coverage_ref } from '../modules/nf-core/samtools/main.nf'
include { samtools_sort as samtools_sort_ref         } from '../modules/nf-core/samtools/main.nf'
include { save_analysis_metadata                     } from '../modules/local/meta/main.nf'
include { seqtk_sample                               } from '../modules/nf-core/seqtk/main.nf'
include { ska_build                                  } from '../modules/nf-core/ska/main.nf'
include { skesa                                      } from '../modules/nf-core/skesa/main.nf'
include { sourmash                                   } from '../modules/nf-core/sourmash/main.nf'
include { spades as spades_illumina                  } from '../modules/nf-core/spades/main.nf'
include { spades as spades_iontorrent                } from '../modules/nf-core/spades/main.nf'

workflow CALL_BACTERIAL_BASE {
    take:
        core_loci_bed
        reference_genome
        reference_genome_dir
        reference_genome_idx
        input_samples
        kraken_db
        target_sample_size
    
    main:
        ch_versions = Channel.empty()

        // Create channel for sample metadata
        Channel.fromPath(input_samples)
            .splitCsv(header:true)
            .tap{ ch_raw_input }
            .map{ row -> get_meta(row) }
            .set{ ch_meta }

        // Create channel for reads
        ch_raw_input
            .map{ row -> get_reads(row) }
            .set{ ch_raw_reads }        

        if ( params.use_hostile ) {
            // remove human reads
            hostile( ch_raw_reads ).reads.set{ ch_depleted_reads }
            ch_versions = ch_versions.mix(hostile.out.versions)
        } else {
            ch_raw_reads.set{ ch_depleted_reads }
        }

        if ( params.target_sample_size ) {
            // downsample reads
            seqtk_sample( ch_depleted_reads, target_sample_size ).reads.set{ ch_depleted_sampled_reads }
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

        Channel.fromPath(params.csv).splitCsv(header:true)
            .map{ row -> get_seqrun_meta(row) }
            .tap{ ch_seqrun_meta }
            .map{ id, sequencing_run, lims_id, sample_name -> [ id, lims_id, sample_name ]}
            .set{ ch_id_meta }

        // create empty channel containing only sample_id
        ch_reads.map{ sample_id, reads -> [ sample_id, [] ] }.set{ ch_empty }

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
        quast(ch_assembly, reference_genome)

        // qc processing - short read
        fastqc(ch_reads_w_meta)
        bwa_mem_ref(ch_reads_w_meta, reference_genome_dir)

        post_align_qc(bwa_mem_ref.out.bam, reference_genome, core_loci_bed)

        // qc processing - long read
        nanoplot(ch_reads_w_meta)

        minimap2_align_ref(ch_reads_w_meta, reference_genome_idx)
        samtools_sort_ref(minimap2_align_ref.out.sam)
        samtools_coverage_ref(samtools_sort_ref.out.bam)

        samtools_sort_ref.out.bam.mix(bwa_mem_ref.out.bam).set{ ch_ref_bam }

        samtools_index_ref(ch_ref_bam)

        if ( params.use_kraken ) {
            kraken(ch_reads, kraken_db)
            bracken(kraken.out.report, kraken_db)
            bracken.out.output.set{ ch_kraken }
            ch_versions = ch_versions.mix(kraken.out.versions)
            ch_versions = ch_versions.mix(bracken.out.versions)
        } else {
            ch_empty.set{ ch_kraken }
        }

        // signature
        sourmash(ch_assembly)
        ska_build(ch_reads)

        ch_versions = ch_versions.mix(bwa_mem_ref.out.versions)
        ch_versions = ch_versions.mix(fastqc.out.versions)
        ch_versions = ch_versions.mix(flye.out.versions)
        ch_versions = ch_versions.mix(medaka.out.versions)
        ch_versions = ch_versions.mix(minimap2_align_ref.out.versions)
        ch_versions = ch_versions.mix(nanoplot.out.versions)
        ch_versions = ch_versions.mix(quast.out.versions)
        ch_versions = ch_versions.mix(samtools_index_ref.out.versions)
        ch_versions = ch_versions.mix(skesa.out.versions)
        ch_versions = ch_versions.mix(ska_build.out.versions)
        ch_versions = ch_versions.mix(sourmash.out.versions)
        ch_versions = ch_versions.mix(spades_illumina.out.versions)
        ch_versions = ch_versions.mix(spades_iontorrent.out.versions)

    emit:
        assembly        = ch_assembly                       // channel: [ val(meta), path(fasta) ]
        bam             = ch_ref_bam                        // channel: [ val(meta), path(bam) ]
        bai             = samtools_index_ref.out.bai        // channel: [ val(meta), path(bai) ]
        empty           = ch_empty                          // channel: [ val(meta) ]
        fastqc          = fastqc.out.output                 // channel: [ val(meta), path(txt) ]
        id_meta         = ch_id_meta                        // channel: [ val(meta), val(meta), val(meta) ]
        kraken          = ch_kraken                         // channel: [ val(meta), path(txt) ]
        metadata        = save_analysis_metadata.out.meta   // channel: [ val(meta), path(json) ]
        qc              = post_align_qc.out.qc              // channel: [ val(meta), path(fasta) ]
        qc_nano_raw     = nanoplot.out.html                 // channel: [ val(meta), path(html) ]
        qc_nano_cov     = samtools_coverage_ref.out.txt     // channel: [ val(meta), path(txt)]
        quast           = quast.out.qc                      // channel: [ val(meta), path(qc) ]
        reads           = ch_reads                          // channel: [ val(meta), path(json) ]
        reads_w_meta    = ch_reads_w_meta                   // channel: [ val(meta), path(meta) ]
        seqplat_meta    = ch_meta                           // channel: [ val(meta), val(str) ]
        ska_build       = ska_build.out.skf                 // channel: [ val(meta), path(skf) ]
        seqrun_meta     = ch_seqrun_meta                    // channel: [ val(meta), val(json), val(json) ]
        sourmash        = sourmash.out.signature            // channel: [ val(meta), path(signature) ]
        versions        = ch_versions                       // channel: [ versions.yml ]
}
