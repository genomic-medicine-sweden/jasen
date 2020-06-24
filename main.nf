#!/usr/bin/env nextflow

samples_ch = Channel.fromFilePairs("${params.fastq_folder}/*{1,2}.fastq.gz")

process krakenBuild {
  publishDir "${params.outdir}/databases", mode: 'copy'
  when: params.kraken2_build

  """
  mkdir -p ${params.outdir}/databases
  UpdateKrakenDatabases.py ${params.outdir}/databases
  """
}

// Could be applied to concat. Uncertain
process fastqc_readqc{
  input:
  file "lane1dir" from samples_ch

  output:
  file "*_fastqc.{zip,html}" into fastqc_results

  """
  touch trams_fastqc.zip
  #fastqc ${lane1dir} --format fastq --threads ${task.cpus}
  #// See manual for: nano, adapters, contaminants
  """
}

forward_ch = Channel.fromPath("${params.fastq_folder}/*1.fastq.gz")
reverse_ch = Channel.fromPath("${params.fastq_folder}/*2.fastq.gz") 

process lane_concatination{
  input:
  file 'forward_concat.fastq.gz' from forward_ch.collectFile() 
  file 'reverse_concat.fastq.gz' from reverse_ch.collectFile()

  output:
  set 'forward_concat.fastq.gz', 'reverse_concat.fastq.gz' into lane_concat_ch

  """
  """
}

process trimmomatic_trimming{
  input:
  file 'sample_master' from lane_concat_ch

  output:
  file "trim_{front_pair, rev_pair, unpair}.fastq.gz" into trimmed_fastq_assembly
  tuple val("trams"), "trim_{front_pair, rev_pair}.fastq.gz" into trimmed_fastq_ref, trimmed_fastq_cont
  
  """
  touch trim_front_pair.fastq.gz trim_front_unpair.fastq.gz trim_rev_pair.fastq.gz trim_rev_unpair.fastq.gz
  cat trim_front_unpair.fastq.gz trim_rev_unpair.fastq.gz >> trim_unpair.fastq.gz
  #trimmomatic PE -threads ${task.cpus} -phred33 $baseDir/${sample_master} $baseDir/ACC6542A4_trim_front_pair.fastq.gz $baseDir/ACC6542A4_trim_front_unpair.fastq.gz $baseDir/ACC6542A4_trim_rev_pair.fastq.gz $baseDir/ACC6542A4_trim_rev_unpair.fastq.gz      ILLUMINACLIP:/home/proj/bin/conda/envs/P_microSALT/share/trimmomatic-0.39-1/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
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
  file(trimmed) from trimmed_fastq_assembly

  output:
  file 'contigs.fasta' into assembled_ch

  """
  touch contigs.fasta
  #spades.py --threads 8 --careful --memory 64 -o /home/proj/production/microbial/results//ACC6417_2020.2.24_11.28.10/ACC6417A161/assembly -1 /home/proj/production/microbial/results//ACC6417_2020.2.24_11.28.10/ACC6417A161/trimmed/ACC6417A161_trim_front_pair.fastq.gz -2 /home/proj/production/microbial/results//ACC6417_2020.2.24_11.28.10/ACC6417A161/trimmed/ACC6417A161_trim_rev_pair.fastq.gz -s /home/proj/production/microbial/results//ACC6417_2020.2.24_11.28.10/ACC6417A161/trimmed/ACC6417A161_trim_unpair.fastq.gz
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
