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
  mkdir -p ${params.krakendb}
  cd ${params.krakendb} && wget ${params.krakendb_url} -O krakendb.tgz
  dlsuf =  tar -tf krakendb.tgz | head -n 1 | tail -c 2
  if [ -f "${params.reference}.sa" ]; then
    tar -xvzf krakendb.tgz --strip 1
  else
    tar -xvzf krakendb.tgz
  fi
  rm krakendb.tgz
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
  tuple "trim_front_pair.fastq.gz", "trim_rev_pair.fastq.gz", "trim_unpair.fastq.gz" into (trimmed_fastq_assembly, trimmed_fastq_ref, trimmed_fastq_cont)
  
  """
  trimmomatic PE -threads ${task.cpus} -phred33 ${forward} ${reverse} trim_front_pair.fastq.gz trim_front_unpair.fastq.gz  trim_rev_pair.fastq.gz trim_rev_unpair.fastq.gz ILLUMINACLIP:${params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  cat trim_front_unpair.fastq.gz trim_rev_unpair.fastq.gz >> trim_unpair.fastq.gz
  """

}
process kraken2_decontamination{
  publishDir "${params.outdir}/kraken2", mode: 'copy', overwrite: true

  input:
    tuple forward, reverse, unpaired from trimmed_fastq_cont

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

process multiqc_report{
  publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: true

  //More inputs as tracks are added
  input:
  file(qreport) from quast_result_ch
  file(freport) from fastqc_results 

  """
  multiqc ${params.outdir} -o ${params.outdir}/multiqc 
  """
}
