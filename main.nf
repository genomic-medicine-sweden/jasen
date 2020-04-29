#!/usr/bin/env nextflow

OUTDIR = params.outdir+'/'+params.subdir

// Get CSV input
//csv = file(params.csv)
Channel
    .fromPath(params.csv).splitCsv(header:true)
    .map{ row-> tuple(row.id, row.species, row.platform, file(row.read1), file(row.read2)) }
    .into { fastq_bwa; fastq_spades; fastq_kraken; fastq_ariba }


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
	sambamba markdup -t ${task.cpus} $bam ${id}.dedup.bam
	"""
}

process kraken {
	publishDir "${OUTDIR}/kraken", mode: 'copy', overwrite: true
	cpus params.cpu_many
	memory '48 GB'
	time '1h'

	input:
		set id, species, platform, file(fastq_r1), file(fastq_r2) from fastq_kraken

	output:
		set id, species, platform, file("${id}.kraken")

	script:
		read_params = fastq_r2.name == 'SINGLE_END' ? 
			"$fastq_r1" :
			"--paired $fastq_r1 $fastq_r2"

	"""
	kraken --fastq-input --check-names --preload \\
		--gzip-compressed \\
		--db ${params.krakendb} \\
		--threads ${task.cpus} \\
		--output kraken.out \\
		$read_params 

	kraken-report --db ${params.krakendb} kraken.out > kraken.rep

	est_abundance.py -k ${params.brackendb} -i kraken.rep -o ${id}.kraken
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
		set id, species, platform, file("${id}.spades.fasta") into asm_quast, asm_mlst, asm_chewbbaca

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


process quast {
	publishDir "${OUTDIR}/qc", mode: 'copy', overwrite: true
	cpus 1
	memory '8 GB'
	time '1h'

	input:
		set id, species, platform, file(asm_fasta) from asm_quast

	output:
		set id, species, platform, file("${id}.quast.tsv")

	script:
		fasta_ref = params.refpath+'/species/'+species+'/ref.fasta'

	"""
	quast.py $asm_fasta -R $fasta_ref -o quast_outdir
 	cp quast_outdir/transposed_report.tsv ${id}.quast.tsv
	"""
}


process mlst {
	publishDir "${OUTDIR}/mlst", mode: 'copy', overwrite: true
	cpus 1
	memory '8 GB'
	time '1h'

	input:
		set id, species, platform, file(asm_fasta) from asm_mlst

	output:
		set id, species, platform, file("${id}.mlst.json"), file("${id}.mlst.novel")


	"""
	mlst --scheme ${species} \\
		--json ${id}.mlst.json --novel ${id}.mlst.novel \\
		${asm_fasta}
	"""
}

process ariba {
	publishDir "${OUTDIR}/ariba", mode: 'copy', overwrite: true
	cpus 1
	memory '8 GB'
	time '1h'

	input:
		set id, species, platform, file(fastq_r1), file(fastq_r2) from fastq_ariba

	output:
		set id, species, platform, file("${id}.ariba.json")

	when:
		fastq_r2 != "SINGLE_END"

	"""
	ariba run --force --threads ${task.cpus} ${params.refpath}/species/${species}/ariba \\
		$fastq_r1 $fastq_r2 ariba.outdir

	ariba summary --col_filter n --row_filter n ariba.summary ariba.outdir/report.tsv

	ariba2json.pl ${params.refpath}/species/${species}/ariba/02.cdhit.all.fa \\
		 ariba.summary.csv ariba.outdir/report.tsv > ${id}.ariba.json
	"""
}

process chewbbaca {
	publishDir "${OUTDIR}/chewbbaca", mode: 'copy', overwrite: true
	cpus 1
	memory '8 GB'
	time '1h'

	input:
		set id, species, platform, file(asm_fasta) from asm_chewbbaca

	output:
		set id, species, platform, file("${id}.chewbbaca")

	script:
		cgmlst_db = params.refpath+'/species/'+species+'/cgmlst'

	"""
	chewBBACA.py AlleleCall -i ${asm_fasta} -g $cgmlst_db --ptf Staphylococcus_aureus.trn --cpu ${task.cpus} -fr \\
	    -o chewbbaca.folder

	sed -e "s/NIPHEM/-/g" -e "s/NIPH/-/g" -e "s/LNF/-/g" -e "s/INF-*//g" -e "s/PLOT[^\t]*/-/g" -e "s/ALM/-/g" -e "s/ASM/-/g" \\
	    chewbbaca.folder/results_alleles.tsv > ${id}.chewbbaca
	"""
}
