#!/usr/bin/env nextflow

process kraken_indexdb {
  publishDir "${params.outdir}/databases", mode: 'copy'
  when: params.kraken2_build

  """
  mkdir -p ${params.outdir}/databases
  UpdateKrakenDatabases.py ${params.outdir}/databases
  """
}

samples_ch = Channel.fromPath("${params.input}/*.fastq.gz")

process fastqc_readqc{
  input:
  file lane1dir from samples_ch

  output:
  file "*_fastqc.{zip,html}" into fastqc_results

  """
  fastqc ${lane1dir} --format fastq --threads ${task.cpus}
  """
}

forward_ch = Channel.fromPath("${params.input}/*1*.fastq.gz")
reverse_ch = Channel.fromPath("${params.input}/*2*.fastq.gz") 

process lane_concatination{
  input:
  file 'forward_concat.fastq.gz' from forward_ch.collectFile() 
  file 'reverse_concat.fastq.gz' from reverse_ch.collectFile()

  output:
  set 'forward_concat.fastq.gz', 'reverse_concat.fastq.gz' into lane_concat_ch

  """
  #Concatination is done via process flow
  """
}

process trimmomatic_trimming{
  input:
  set forward, reverse from lane_concat_ch

  output:
  file "trim_{front_pair, rev_pair, unpair}.fastq.gz" into trimmed_fastq_assembly
  tuple val("trams"), "trim_{front_pair, rev_pair}.fastq.gz" into trimmed_fastq_ref, trimmed_fastq_cont
  
  """
  trimmomatic PE -threads ${task.cpus} -phred33 ${forward} ${reverse} trim_front_pair.fastq.gz trim_rev_pair.fastq.gz trim_front_unpair.fastq.gz trim_rev_unpair.fastq.gz ILLUMINACLIP:${baseDir}/assets/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  cat trim_front_unpair.fastq.gz trim_rev_unpair.fastq.gz >> trim_unpair.fastq.gz
  """

}

process kraken2_decontamination{
  publishDir "${params.outdir}/kraken2", mode: 'copy'

  input:
    set val(name), file(reads) from trimmed_fastq_cont

  output:
    set val(name), file("*.kraken.out") into kraken_out
    set val(name), file("*.kraken.report") into kraken_report 


  script:
    out = name+".kraken.out"
    kreport = name+".kraken.report"
    if (params.pairedEnd){
      """
      touch ${out} ${kreport}
      #kraken2 --db ${params.krakendb} --threads ${task.cpus} --output $out --report $kreport --paired ${reads[0]} ${reads[1]}
      """    
    } else {
      """
      touch ${out} ${kreport}
      #kraken2 --db ${params.krakendb} --threads ${task.cpus} --output $out --report $kreport ${reads[0]}
      """
    }
}


process spades_assembly{
  input:
  set val(name), file(trimmed) from trimmed_fastq_ass

  output:
  file 'spades/contigs.fasta' into assembled_ch

  script:
  """
  spades.py --threads ${task.cpus} --careful -o spades -1 trimmed[0] -2 trimmed[1] -s trimmed[2]
  """
}

process quast_assembly_qc{
  input:
  file contig from assembled_ch 

  """
  #quast.py /home/proj/production/microbial/results//ACC6417_2020.2.24_11.28.10/ACC6417A161/assembly/contigs.fasta -o /home/proj/production/microbial/results//ACC6417_2020.2.24_11.28.10/ACC6417A161/assembly/quast
 """
}

process bwa_read_mapping{
  input:
  file(trimmed) from trimmed_fastq_ref

  """
  #bwa mem -M -t 8 /home/proj/production/microbial/references/genomes/NC_000913.fasta /home/proj/production/microbial/results//ACC6417_2020.2.24_11.28.10/ACC6417A161/trimmed/ACC6417A161_trim_front_pair.fastq.gz /home/proj/production/microbial/results//ACC6417_2020.2.24_11.28.10/ACC6417A161/trimmed/ACC6417A161_trim_rev_pair.fastq.gz > /home/proj/production/microbial/results//ACC6417_2020.2.24_11.28.10/ACC6417A161/alignment/ACC6417A161_NC_000913.sam
  """
}
