#!/usr/bin/env nextflow

process bwa_index_reference{
  cpus 1

  output:
  file "bypass" into bwa_indexes_ch

  """
  if [ ! -f "${params.reference}.sa" ]; then
    bwa index ${params.reference}
  fi
  touch bypass
  """
}

process kraken2_db_download{
  cpus 1

  output:
  file 'database.rdy' into kraken_init_ch

  """
  if ${params.kraken_db_download} ; then
    wd=\$(pwd)
    export PATH=$PATH:$baseDir/bin/
    mkdir -p ${params.krakendb}
    cd ${params.krakendb} && wget ${params.krakendb_url} -O krakendb.tgz
    dlsuf=`tar -tf krakendb.tgz | head -n 1 | tail -c 2`
    if [ -f "${params.reference}.sa" ]; then
      tar -xvzf krakendb.tgz --strip 1
    else
      tar -xvzf krakendb.tgz
    fi
    rm krakendb.tgz
    cd \${wd} && touch database.rdy
  else
    cd \${wd} && touch database.rdy
  fi
  
  """
}

process ariba_db_download{

  output:
  file 'database.rdy' into ariba_init_ch

  """
  if  ${params.ariba_db_download} ; then
    ariba getref resfinder resfinder
    ariba prepareref --force -f ./resfinder.fa -m ./resfinder.tsv --threads ${task.cpus} ${params.aribadb}
    mv resfinder.fa ${params.aribadb}
    mv resfinder.tsv ${params.aribadb}
    touch database.rdy
  else
    touch database.rdy
  fi
  """

}

samples_ch = Channel.fromPath("${params.input}/*.{fastq.gz,fsa.gz,fa.gz,fastq,fsa,fa}")

process fastqc_readqc{
  publishDir "${params.outdir}/fastqc", mode: 'copy', overwrite: true

  input:
  file lane1dir from samples_ch

  output:
  file "*_fastqc.html" into fastqc_results

  """
  fastqc ${params.input}/${lane1dir} --format fastq --threads ${task.cpus} -o .
  """
}

forward_ch = Channel.fromPath("${params.input}/*1*.{fastq.gz,fsa.gz,fa.gz,fastq,fsa,fa}")
reverse_ch = Channel.fromPath("${params.input}/*2*.{fastq.gz,fsa.gz,fa.gz,fastq,fsa,fa}")
 

process lane_concatination{
  publishDir "${params.outdir}/concatinated", mode: 'copy', overwrite: true
  cpus 1

  input:
  file 'forward_concat.fastq.gz' from forward_ch.collectFile() 
  file 'reverse_concat.fastq.gz' from reverse_ch.collectFile()

  output:
  tuple 'forward_concat.fastq.gz', 'reverse_concat.fastq.gz' into lane_concat_ch

  """
  #Concatination is done via process flow
  """
}

process trimmomatic_trimming{
  publishDir "${params.outdir}/trimmomatic", mode: 'copy', overwrite: true

  input:
  tuple forward, reverse from lane_concat_ch

  output:
  tuple "trim_front_pair.fastq.gz", "trim_rev_pair.fastq.gz", "trim_unpair.fastq.gz" into (trimmed_fastq_assembly, trimmed_fastq_ref, trimmed_fastq_cont, trimmed_ariba_cont)
  
  """
  trimmomatic PE -threads ${task.cpus} -phred33 ${forward} ${reverse} trim_front_pair.fastq.gz trim_front_unpair.fastq.gz  trim_rev_pair.fastq.gz trim_rev_unpair.fastq.gz ILLUMINACLIP:${params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  cat trim_front_unpair.fastq.gz trim_rev_unpair.fastq.gz >> trim_unpair.fastq.gz
  """

}

process ariba_resistancefind{
  publishDir "${params.outdir}/ariba", mode: 'copy', overwrite: true

  input:
  tuple forward, reverse, unpaired from trimmed_ariba_cont 
  file(database_initalization) from ariba_init_ch 

  output:
  file 'resistance.summary' into ariba_output
  

  """
  ariba run --spades_options careful --force --threads ${task.cpus} ${params.aribadb} ${forward} ${reverse} ${params.outdir}/ariba
  """
}

process ariba_stats{
  publishDir "${params.outdir}/ariba", mode: 'copy', overwrite: true
  cpus 1

  input:
  file(aribasummary) from ariba_output

  output:
  file 'report.tsv' into ariba_summary_output 

  """
  ariba summary --col_filter n --row_filter n resistance.summary ${params.outdir}/ariba/report.tsv
  """
}

process kraken2_decontamination{
  publishDir "${params.outdir}/kraken2", mode: 'copy', overwrite: true

  input:
  tuple forward, reverse, unpaired from trimmed_fastq_cont
  file(db_initialized) from kraken_init_ch


  output:
  file "kraken.out" into kraken_out
  file "kraken.report" into kraken_report 


  """
  kraken2 --db ${params.krakendb} --threads ${task.cpus} --output kraken.out --report kraken.report --paired ${forward} ${reverse}
  """    
}
process spades_assembly{
  publishDir "${params.outdir}/spades", mode: 'copy', overwrite: true

  input:
  file(reads) from trimmed_fastq_assembly

  output:
  file 'contigs.fasta' into (assembled_ch, mlst_ch)

  script:
  """
  spades.py --threads ${task.cpus} --careful -o . -1 ${reads[0]} -2 ${reads[1]} -s ${reads[2]}
  """
}

process mlst_lookup{
  publishDir "${params.outdir}/mlst", mode: 'copy', overwrite: true

  input:
  file contig from mlst_ch


  """
  mlst $contig --threads ${task.cpus} --json mlst.json --novel novel.fasta --minid 99.5 --mincov 95
  """
}

process quast_assembly_qc{
  publishDir "${params.outdir}/quast", mode: 'copy', overwrite: true
  cpus 1

  input:
  file contig from assembled_ch 

  output:
  file 'report.tsv' into quast_result_ch

  """
  quast.py $contig -o .
 """
}

process bwa_read_mapping{
  publishDir "${params.outdir}/bwa", mode: 'copy', overwrite: true

  input:
  file(trimmed) from trimmed_fastq_ref
  file(bypass) from bwa_indexes_ch

  output:
  file 'alignment.sam' into bwa_mapping_ch

  """
  bwa mem -M -t ${task.cpus} ${params.reference} ${trimmed[0]} ${trimmed[1]} > alignment.sam
  """
}

process samtools_bam_conversion{
  publishDir "${params.outdir}/bwa", mode: 'copy', overwrite: true

  input:
  file(aligned_sam) from bwa_mapping_ch

  output:
  file 'alignment_sorted.bam' into bwa_sorted_ch, bwa_sorted_ch2

  """
  samtools view --threads ${task.cpus} -b -o alignment.bam -T ${params.reference} ${aligned_sam}
  samtools sort --threads ${task.cpus} -o alignment_sorted.bam alignment.bam
  

  """
}

process samtools_duplicates_stats{
  publishDir "${params.outdir}/samtools", mode: 'copy', overwrite: true

  input:
  file(align_sorted) from bwa_sorted_ch2

  output:
  tuple 'samtools_map.txt', 'samtools_raw.txt' into samtools_dup_results

  """
  samtools flagstat ${align_sorted} &> samtools_map.txt
  samtools view -c ${align_sorted} &> samtools_raw.txt
  """
}

process picard_markduplicates{
  publishDir "${params.outdir}/picard", mode: 'copy', overwrite: true
  cpus 1

  input:
  file(align_sorted) from bwa_sorted_ch

  output:
  file 'alignment_sorted_rmdup.bam' into bam_rmdup_ch, bam_rmdup_ch2
  file 'picard_dupstats.txt' into picard_stats_hist_ch

  """
  picard MarkDuplicates I=${align_sorted} O=alignment_sorted_rmdup.bam M=picard_dupstats.txt REMOVE_DUPLICATES=true
  """
}



process picard_qcstats{
  publishDir "${params.outdir}/picard", mode: 'copy', overwrite: true
  cpus 1

  input:
  file(alignment_sorted_rmdup) from bam_rmdup_ch
  
  output:
  tuple 'picard_stats.txt', 'picard_insstats.txt' into picard_stats_ch

  """
  picard CollectInsertSizeMetrics I=${alignment_sorted_rmdup} O=picard_stats.txt H=picard_insstats.txt

  """
}

process samtools_deduplicated_stats{
  publishDir "${params.outdir}/samtools", mode: 'copy', overwrite: true

  input:  
  file(alignment_sorted_rmdup) from bam_rmdup_ch2

  output:
  tuple 'samtools_ref.txt', 'samtools_cov.txt' into samtools_dedup_results

  """
  samtools index ${alignment_sorted_rmdup}
  samtools idxstats ${alignment_sorted_rmdup} &> samtools_ref.txt
  samtools stats --coverage 1,10000,1 ${alignment_sorted_rmdup} |grep ^COV | cut -f 2- &> samtools_cov.txt

  """

}

process multiqc_report{
  publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: true
  cpus 1

  //More inputs as tracks are added
  input:
  file(qreport) from quast_result_ch
  file(freport) from fastqc_results 

  """
  multiqc ${params.outdir} -o ${params.outdir}/multiqc 
  """
}
