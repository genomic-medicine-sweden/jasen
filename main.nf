#!/usr/bin/env nextflow

OUTDIR = params.outdir+'/'+params.subdir

// Get CSV input
//csv = file(params.csv)
Channel
    .fromPath(params.csv).splitCsv(header:true)
    .map{ row-> tuple(row.id, row.species, row.platform, file(row.read1), file(row.read2)) }
    .into { fastq_bwa; fastq_spades }



process bwa_align {
	cpus params.cpu_bwa
	memory '32 GB'
	time '1h'
	    
	input: 
		set id, species, platform, file(fastq_r1), file(fastq_r2) from fastq_bwa

	output:
		set id, species, platform, file("${id}.bwa.sort.bam"), file("${id}.bwa.sort.bam.bai") into bam_markdup

	script:
		fasta_ref = params.refpath+'/species/'+species+'/ref.fasta'
		read2 = fastq_r2.name == 'SINGLE_END' ? '' : "$fastq_r2"

		"""
		bwa mem -R '@RG\\tID:${id}\\tSM:${id}\\tPL:illumina' -M -t ${task.cpus} $fasta_ref $fastq_r1 $read2 \\
		| samtools view -Sb - \\
		| samtools sort -o ${id}.bwa.sort.bam -

		samtools index ${id}.bwa.sort.bam
		"""
}


process bam_markdup {
	publishDir "${OUTDIR}/bam", mode: 'copy', overwrite: true
	cpus params.cpu_many
	memory '32 GB'
	time '1h'
    
	input:
		set id, species, platform, file(bam), file(bai) from bam_markdup

	output:
		set id, species, platform, file("${id}.dedup.bam"), file("${id}.dedup.bam.bai") into bam_qc

	"""
	sambamba markdup -t ${task.cpu} $bam ${id}.dedup.bam
	"""
}


process spades_assembly {
	publishDir "${OUTDIR}/assembly", mode: 'copy', overwrite: true
	cpus params.cpu_spades
	memory '16 GB'
	time '2h'

	input:
		set id, species, platform, file(fastq_r1), file(fastq_r2) from fastq_spades

	output:
		set id, species, platform, file("${id}.spades.fasta")

	script:

		opt_platform = platform == 'iontorrent' ? '--iontorrent --careful' : '--only-assembler'
		opt_reads = fastq_r2.name != 'SINGLE_END' ? "-1 $fastq_r1 -2 $fastq_r2" : "-s $fastq_r1"

		"""
		spades.py -k 21,33,55,77 -t ${task.cpus} \\
			$opt_platform \\
			$opt_reads \\
			-o ${id}.spades.fasta
		"""
}


