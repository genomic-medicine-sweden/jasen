#!/usr/bin/env nextflow

process bwa_index_reference{

  output:
  file "bypass" into bwa_indexes_ch

  """
  if [ ! -f "${params.reference}.sa" ]; then
    bwa index ${params.reference}
  fi
  touch bypass
  """
}

process kraken_db_download {
  when: params.kraken_db_download

  """
  export PATH=$PATH:$baseDir/bin/
  kpath=dirname ${params.krakendb}
  mkdir -p ${kpath}
  cd ${kpath} && wget ${krakendb_url}
  """
}


samples_ch = Channel.fromPath("${params.input}/*.fastq.gz")

process fastqc_readqc{
  input:
  file lane1dir from samples_ch

  output:
  file "*_fastqc.{zip,html}" into fastqc_results

  """
  fastqc ${params.input}/${lane1dir} --format fastq --threads ${task.cpus} -o .
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
  tuple "trim_front_pair.fastq.gz", "trim_rev_pair.fastq.gz", "trim_unpair.fastq.gz" into (trimmed_fastq_assembly, trimmed_fastq_ref)
  tuple val("trams"), "trim_front_pair.fastq.gz", "trim_rev_pair.fastq.gz" into trimmed_fastq_cont
  
  """
  trimmomatic PE -threads ${task.cpus} -phred33 ${forward} ${reverse} trim_front_pair.fastq.gz trim_front_unpair.fastq.gz  trim_rev_pair.fastq.gz trim_rev_unpair.fastq.gz ILLUMINACLIP:${params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  cat trim_front_unpair.fastq.gz trim_rev_unpair.fastq.gz >> trim_unpair.fastq.gz
  """

}
process kraken2_decontamination{
  input:
    set val(name), file(reads) from trimmed_fastq_cont

  output:
    set val(name), file("kraken.out") into kraken_out
    set val(name), file("kraken.report") into kraken_report 


  """
  kraken2 --db ${params.krakendb} --threads ${task.cpus} --output kraken.out --report kraken.report --paired ${reads[0]} ${reads[1]}
  """    
}
process spades_assembly{
  input:
  file(reads) from trimmed_fastq_assembly

  output:
  file 'spades/contigs.fasta' into (assembled_ch, mlst_ch)

  script:
  """
  spades.py --threads ${task.cpus} --careful -o spades -1 ${reads[0]} -2 ${reads[1]} -s ${reads[2]}
  """
}

process mlst_lookup{
  input:
  file contig from mlst_ch


  """
  mlst $contig --threads ${task.cpus} --json mlst.json --novel novel.fasta --minid 99.5 --mincov 95
  """
}

process quast_assembly_qc{
  input:
  file contig from assembled_ch 

  output:
  file 'report.tsv' into quast_result_ch

  """
  quast.py $contig -o .
 """
}

process bwa_read_mapping{
  input:
  file(trimmed) from trimmed_fastq_ref
  file(bypass) from bwa_indexes_ch

  output:
  file 'alignment.sam' into bwa_mapping_ch

  """
  bwa mem -M -t ${task.cpus} ${params.reference} ${trimmed[0]} ${trimmed[1]} > alignment.sam
  """
}
