#!/usr/bin/env nextflow


if (!(params.pkm && params.location)) {
  exit 1, "YOU HAVE TO PROVIDE A LOCATION AND PACKAGE MANAGER PROFILE E.g. 'nextflow run main.nf -profile local,conda'"
}

process setup_workdirectories{
  label 'min_allocation'

  output:
  file "assets.rdy" into assets_done

  """
  cp -r ${baseDir}/assets/* ${params.assets}
  mkdir -p ${params.rootdir} 
  mkdir -p ${params.work}
  touch assets.rdy
  """
}

process bwa_index_reference{
  label 'min_allocation'

  input:
  file assets_rdy from assets_done

  output:
  file "database.rdy" into bwa_indexes

  """
  if [ ! -f "${params.reference}.sa" ]; then
    bwa index ${params.reference}
  fi
  touch database.rdy
  """
}

process cgmlst_db_init{
  label 'min_allocation'

  output:
  file 'database.rdy' into chewie_init
  file 'chewiedb.zip' into chewie_source

  """
  if ${params.chewbbaca_db_download} ; then
    export PATH=\$PATH:$baseDir/bin/
    mkdir -p ${params.chewbbacadb} 
    wget ${params.chewbbacadb_url} -O chewiedb.zip
    unzip -o chewiedb.zip -d ${params.chewbbacadb} 
    chewBBACA.py PrepExternalSchema -i ${params.chewbbacadb} -o ${params.chewbbacadb}/schema --cpu ${task.cpus}
    touch database.rdy
  else
    touch database.rdy
  fi
  """
}

process kraken2_db_download{
  label 'min_allocation'

  output:
  file 'database.rdy' into kraken2_init
  file 'krakendb.tgz' into kraken2_source

  """
  if ${params.kraken_db_download} ; then
    export PATH=\$PATH:$baseDir/bin/
    mkdir -p ${params.krakendb}
    wget ${params.krakendb_url} -O krakendb.tgz
    # dlsuf=`tar -tf krakendb.tgz | head -n 1 | tail -c 2`
    if [ -f "${params.reference}.sa" ]; then
      tar -xvzf krakendb.tgz -C ${params.krakendb} --strip 1
    else
      tar -xvzf krakendb.tgz -C ${params.krakendb}
    fi
    # rm krakendb.tgz
    touch database.rdy
  else
    touch database.rdy
  fi
  """
}

process ariba_db_download{
  label 'modest_allocation'
 
  output:
  file 'database.rdy' into ariba_init

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

samples = Channel.fromPath("${params.input}/*.{fastq.gz,fsa.gz,fa.gz,fastq,fsa,fa}")

process fastqc_readqc{
  label 'modest_allocation'

  publishDir "${params.outdir}/fastqc", mode: 'copy', overwrite: true

  input:
  file lane1dir from samples

  output:
  file "*_fastqc.html" into fastqc_results

  """
  fastqc ${params.input}/${lane1dir} --format fastq --threads ${task.cpus} -o .
  """
}

forward = Channel.fromPath("${params.input}/*_R1_*.{fastq.gz,fsa.gz,fa.gz,fastq,fsa,fa}")
reverse = Channel.fromPath("${params.input}/*_R2_*.{fastq.gz,fsa.gz,fa.gz,fastq,fsa,fa}")


process lane_concatination{
  label 'min_allocation'

  publishDir "${params.outdir}/concatinated", mode: 'copy', overwrite: true

  input:
  file 'forward_concat.fastq.gz' from forward.collectFile()
  file 'reverse_concat.fastq.gz' from reverse.collectFile()

  output:
  tuple 'forward_concat.fastq.gz', 'reverse_concat.fastq.gz' into lane_concat

  """
  #Concatination is done via process flow
  """
}

process trimmomatic_trimming{
  label 'min_allocation'

  publishDir "${params.outdir}/trimmomatic", mode: 'copy', overwrite: true

  input:
  tuple forward, reverse from lane_concat

  output:
  tuple "trim_front_pair.fastq.gz", "trim_rev_pair.fastq.gz", "trim_unpair.fastq.gz" into (trimmed_sample_1, trimmed_sample_2, trimmed_sample_3, trimmed_sample_4)

  """
  trimmomatic PE -threads ${task.cpus} -phred33 ${forward} ${reverse} trim_front_pair.fastq.gz trim_front_unpair.fastq.gz  trim_rev_pair.fastq.gz trim_rev_unpair.fastq.gz ILLUMINACLIP:${params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  cat trim_front_unpair.fastq.gz trim_rev_unpair.fastq.gz >> trim_unpair.fastq.gz
  """

}

process ariba_resistancefind{
  label 'modest_allocation'

  publishDir "${params.outdir}/ariba", mode: 'copy', overwrite: true, pattern: 'motif_report.tsv'

  input:
  tuple forward, reverse, unpaired from trimmed_sample_4
  file(database_initalization) from ariba_init

  output:
  file 'motif_report.tsv' into ariba_output


  """
  ariba run --spades_options careful --force --threads ${task.cpus} ${params.aribadb} ${forward} ${reverse} outdir
  mv outdir/report.tsv motif_report.tsv
  """
}

process ariba_stats{
  label 'min_allocation'

  publishDir "${params.outdir}/ariba", mode: 'copy', overwrite: true
  cpus 1

  input:
  file(report) from ariba_output

  output:
  tuple 'summary.csv', 'motif_report.json' into ariba_summary_output

  """
  ariba summary --col_filter n --row_filter n summary ${report} 
  python3 $baseDir/bin/tsv_to_json.py ${report} motif_report.json 
  """
}

process kraken2_decontamination{
  label 'max_allocation'

  publishDir "${params.outdir}/kraken2", mode: 'copy', overwrite: true

  input:
  tuple forward, reverse, unpaired from trimmed_sample_3
  file(db_initialized) from kraken2_init


  output:
  tuple "kraken_out.tsv", "kraken_report.tsv" into kraken2_output


  """
  kraken2 --db ${params.krakendb} --threads ${task.cpus} --output kraken_out.tsv --report kraken_report.tsv --paired ${forward} ${reverse}
  """
}
process spades_assembly{
  label 'max_allocation'

  publishDir "${params.outdir}/spades", mode: 'copy', overwrite: true

  input:
  file(reads) from trimmed_sample_1

  output:
  file 'scaffolds.fasta' into (assembled_sample_1, assembled_sample_2, assembled_sample_3)

  script:
  """
  spades.py --threads ${task.cpus} --careful -o . -1 ${reads[0]} -2 ${reads[1]} -s ${reads[2]}
  """
}

process mlst_lookup{
  label 'min_allocation'

  publishDir "${params.outdir}/mlst", mode: 'copy', overwrite: true

  input:
  file contig from assembled_sample_1

  output:
  file 'mlst.json' into mlst_output

  """
  mlst $contig --threads ${task.cpus} --json mlst.json --novel novel_mlst.fasta --minid 99.5 --mincov 95
  """
}

process chewbbaca_cgmlst{
  label 'max_allocation'
  publishDir "${params.outdir}/cgmlst", mode: 'copy', overwrite: true

  input:
  file contig from assembled_sample_3
  file 'database.dry' from chewie_init

  output:
  tuple 'cgmlst_alleles.json', 'cgmlst_stats.json' into cgmlst_results

  """
  yes | chewBBACA.py AlleleCall --fr -i \${PWD} -g ${params.chewbbacadb}/schema --json --cpu ${task.cpus} -o \${PWD} --ptf ${params.prodigal_file}
  mv results_*/* .
  mv results_alleles.json cgmlst_alleles.json
  mv results_statistics.json cgmlst_stats.json
  """
}


process quast_assembly_qc{
  label 'min_allocation'

  publishDir "${params.outdir}/quast", mode: 'copy', overwrite: true

  input:
  file contig from assembled_sample_2

  output:
  file 'quast_report.tsv' into quast_result, quast_result_2

  """
  quast.py $contig -o . -r ${params.reference} -t ${task.cpus}
  cp report.tsv quast_report.tsv
  """
}

process quast_json_conversion{
  label 'min_allocation'  

  publishDir "${params.outdir}/quast", mode: 'copy', overwrite: true
  cpus 1

  input:
  file(quastreport) from quast_result_2

  output:
  file 'quast_report.json' into quast_result_json

  """
  python3 $baseDir/bin/quast_to_json.py $quastreport quast_report.json
  """
}


process bwa_read_mapping{
  label 'max_allocation'

  publishDir "${params.outdir}/bwa", mode: 'copy', overwrite: true

  input:
  file(trimmed) from trimmed_sample_2
  file(database_initalization) from bwa_indexes

  output:
  file 'alignment.sam' into mapped_sample

  """
  bwa mem -M -t ${task.cpus} ${params.reference} ${trimmed[0]} ${trimmed[1]} > alignment.sam
  """
}

process samtools_bam_conversion{
  label 'min_allocation'

  publishDir "${params.outdir}/bwa", mode: 'copy', overwrite: true

  input:
  file(aligned_sam) from mapped_sample

  output:
  file 'alignment_sorted.bam' into sorted_sample_1, sorted_sample_2

  """
  samtools view --threads ${task.cpus} -b -o alignment.bam -T ${params.reference} ${aligned_sam}
  samtools sort --threads ${task.cpus} -o alignment_sorted.bam alignment.bam
  """
}

process samtools_duplicates_stats{
  label 'min_allocation'

  publishDir "${params.outdir}/samtools", mode: 'copy', overwrite: true

  input:
  file(align_sorted) from sorted_sample_1

  output:
  tuple 'samtools_flagstats.txt', 'samtools_total_reads.txt' into samtools_duplicated_results

  """
  samtools flagstat ${align_sorted} &> samtools_flagstats.txt
  samtools view -c ${align_sorted} &> samtools_total_reads.txt
  """
}

process picard_markduplicates{
  label 'min_allocation'

  publishDir "${params.outdir}/picard", mode: 'copy', overwrite: true
  cpus 1

  input:
  file(align_sorted) from sorted_sample_2

  output:
  file 'alignment_sorted_rmdup.bam' into deduplicated_sample, deduplicated_sample_2, deduplicated_sample_3
  file 'picard_duplication_stats.txt' into picard_histogram_output

  """
  picard MarkDuplicates I=${align_sorted} O=alignment_sorted_rmdup.bam M=picard_duplication_stats.txt REMOVE_DUPLICATES=true
  """
}

process samtools_calling{
  label 'min_allocation'

  publishDir "${params.outdir}/snpcalling", mode: 'copy', overwrite: true

  input:
  file(align_sorted_rmdup) from deduplicated_sample

  output:
  file 'samtools_calls.bam' into called_sample

  """
  samtools view -@ ${task.cpus} -h -q 1 -F 4 -F 256 ${align_sorted_rmdup} | grep -v XA:Z | grep -v SA:Z| samtools view -b - > samtools_calls.bam
  """
}


process vcftools_snpcalling{
  label 'min_allocation'

  publishDir "${params.outdir}/snpcalling", mode: 'copy', overwrite: true

  input:
  file(samhits) from called_sample

  output:
  file 'vcftools.recode.bcf' into snpcalling_output

  """
  vcffilter="--minQ 30 --thin 50 --minDP 3 --min-meanDP 20"
  bcffilter="GL[0]<-500 & GL[1]=0 & QR/RO>30 & QA/AO>30 & QUAL>5000 & ODDS>1100 & GQ>140 & DP>100 & MQM>59 & SAP<15 & PAIRED>0.9 & EPP>3"


  freebayes -= --pvar 0.7 -j -J --standard-filters -C 6 --min-coverage 30 --ploidy 1 -f ${params.reference} -b ${samhits} -v freebayes.vcf
  bcftools view freebayes.vcf -o unfiltered_bcftools.bcf.gz -O b --exclude-uncalled --types snps
  bcftools index unfiltered_bcftools.bcf.gz
  bcftools view unfiltered_bcftools.bcf.gz -i \${bcffilter} -o bcftools.bcf.gz -O b
  vcftools --bcf bcftools.bcf.gz \${vcffilter} --remove-filtered-all --recode-INFO-all --recode-bcf --out vcftools

  """
}

process snp_translation{
  publishDir "${params.outdir}/snpcalling", mode: 'copy', overwrite: true

  label 'min_allocation'

  input:
  file bcf_file from snpcalling_output

  output:
  tuple 'vcftools.recode.vcf', 'snp_report.tsv' into snp_translated_output
  file 'snp_report.json' into snp_json_output

  script:
  """
  bcftools view ${bcf_file} > vcftools.recode.vcf
  gatk VariantsToTable -V vcftools.recode.vcf -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F DP -F I16 -F QS -F MQ0F -GF PL -O snp_report.tsv
  python3 $baseDir/bin/tsv_to_json.py snp_report.tsv snp_report.json
  """

}


process picard_qcstats{
  label 'min_allocation'

  publishDir "${params.outdir}/picard", mode: 'copy', overwrite: true

  input:
  file(alignment_sorted_rmdup) from deduplicated_sample_2

  output:
  tuple 'picard_stats.txt', 'picard_insert_distribution.pdf' into picard_output

  """
  picard CollectInsertSizeMetrics I=${alignment_sorted_rmdup} O=picard_stats.txt H=picard_insert_distribution.pdf

  """
}

process samtools_deduplicated_stats{
  label 'min_allocation'

  publishDir "${params.outdir}/samtools", mode: 'copy', overwrite: true

  input:
  file(alignment_sorted_rmdup) from deduplicated_sample_3

  output:
  tuple 'samtools_idxstats.tsv', 'samtools_coverage_distribution.tsv' into samtools_deduplicated_output

  """
  samtools index ${alignment_sorted_rmdup}
  samtools idxstats ${alignment_sorted_rmdup} &> samtools_idxstats.tsv
  samtools stats --coverage 1,10000,1 ${alignment_sorted_rmdup} |grep ^COV | cut -f 2- &> samtools_coverage_distribution.tsv

  """

}

/*
The following reports are generated ( * = Supported in multiqc):

Ariba summary
* Kraken report
MLST report, MLST novel
* Picard Insert Size
* Picard MarkDuplicates
* Quast report
Quast json
* Samtools flagstat
* Samtools idxstats
Samtools coverage distribution
Samtools total reads
SNPcalling

*/

process multiqc_report{
  label 'min_allocation'

  publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: true

  //More inputs as tracks are added
  input:
  file(quast_report) from quast_result
  file(fastqc_report) from fastqc_results
  tuple snp_vcf, snp_tsv from snp_translated_output
  tuple picard_stats, picard_insert_stats from picard_output
  tuple kraken_output, kraken_report from kraken2_output
  tuple samtools_map, samtools_raw from samtools_duplicated_results

  output:
  file 'multiqc_report.html' into multiqc_output
  file 'multiqc_data/multiqc_data.json' into multiqc_json
  // MultiQC_data contains a lot delimited files. May be useful later

  """
  multiqc ${params.outdir} -f -k json -o \$(pwd)
  """
}

process output_collection{
  label 'min_allocation'
"""
"""
}

process json_collection{
  label 'min_allocation'

  publishDir "${params.outdir}/jsoncollection", mode: 'copy', overwrite: true

  input:
  file (mlstjson) from mlst_output
  file (multiqcjson) from multiqc_json
  file (aribajson) from ariba_summary_output
  file (quastjson) from quast_result_json
  file (snpreport) from snp_json_output
  tuple (cgmlst_res, cgmlst_stats) from cgmlst_results
  
  output:
  tuple 'merged_report.json', mlstjson, multiqcjson, aribajson, quastjson, snpreport, cgmlst_res into json_collection

  """
  touch merged_reports.json
  cat ${mlstjson} >> merged_report.json
  cat ${aribajson} >> merged_report.json
  cat ${quastjson} >> merged_report.json
  cat ${snpreport} >> merged_report.json
  cat ${multiqcjson} >> merged_report.json
  cat ${cgmlst_res} >> merged_report.json
  """
}


