#!/usr/bin/env nextflow

OUTDIR = params.outdir+'/'+params.subdir

// Get CSV input
//csv = file(params.csv)
Channel
    .fromPath(params.csv).splitCsv(header:true)
    .map{ row-> tuple(row.id, row.species, row.platform, file(row.read1), file(row.read2)) }
    .into { fastq_bwa; fastq_spades; fastq_kraken; fastq_ariba; fastq_maskpolymorph; fastq_register }


process bwa_align {
	cpus params.cpu_bwa
	memory '32 GB'
	time '1h'

	input:
		set id, species, platform, file(fastq_r1), file(fastq_r2) from fastq_bwa

	output:
		set id, species, platform, file("${id}.bwa.sort.bam"), file("${id}.bwa.sort.bam.bai") into bam_markdup, bam_postalignqc

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
		set id, species, platform, file("${id}.kraken") into kraken_export

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
		set id, species, platform, file("${id}.fasta") into asm_quast, asm_mlst, asm_chewbbaca, asm_maskpolymorph

	script:

		opt_platform = platform == 'iontorrent' ? '--iontorrent --careful' : '--only-assembler'
		opt_reads = fastq_r2.name != 'SINGLE_END' ? "-1 $fastq_r1 -2 $fastq_r2" : "-s $fastq_r1"

	"""
	spades.py -k 21,33,55,77 -t ${task.cpus} \\
		$opt_platform \\
		$opt_reads \\
		-o spades
	mv spades/contigs.fasta ${id}.fasta
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
		set id, species, platform, file("${id}.quast.tsv") into quast_register, quast_export

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
		set id, species, platform, file("${id}.mlst.json") into mlst_export
		file("${id}.mlst.novel") optional true


	"""
	mlst --scheme ${species} \\
		--json ${id}.mlst.json --novel ${id}.mlst.novel \\
		${asm_fasta}
	"""
}


process ariba {
	publishDir "${OUTDIR}/ariba", mode: 'copy', overwrite: true
	// cache 'deep'
	cpus params.cpu_many
	memory '16 GB'
	time '1h'

	input:
		set id, species, platform, file(fastq_r1), file(fastq_r2) from fastq_ariba

	output:
		set id, species, platform, file("${id}.ariba.json") into ariba_export

	when:
		fastq_r2 != "SINGLE_END"

	"""
	ariba run --force --threads ${task.cpus} ${params.refpath}/species/${species}/ariba \\
		$fastq_r1 $fastq_r2 ariba.outdir

	ariba summary --col_filter n --row_filter n ariba.summary ariba.outdir/report.tsv

	ariba2json.pl ${params.refpath}/species/${species}/ariba/02.cdhit.all.fa \\
		 ariba.summary.csv ariba.outdir/report.tsv > ${id}.ariba.json
	#echo foo > ${id}.ariba.json
	"""
}


process maskpolymorph {
	publishDir "${OUTDIR}/maskepolymorph", mode: 'copy', overwrite: true
	cpus params.cpu_bwa
	memory '32 GB'
	time '1h'
	// cache 'deep'

	input:
		set id, species, platform, file(asm_fasta), file(fastq_r1), file(fastq_r2) from asm_maskpolymorph.join(fastq_maskpolymorph, by:[0,1,2])

	output:
		set id, species, platform, file("${id}.spades.masked") into maskpoly_chewbbaca

	script:
		read2 = fastq_r2.name == 'SINGLE_END' ? '' : "$fastq_r2"

	"""
	bwa index $asm_fasta
	bwa mem -t ${task.cpus} -R '@RG\\tID:${id}\\tSM:${id}\\tPL:${platform}' $asm_fasta $fastq_r1 $read2 \\
		| samtools view -b - -@ ${task.cpus} \\
		| samtools sort -@ ${task.cpus} - -o contigs.fasta.sort.bam
	samtools index contigs.fasta.sort.bam

	freebayes -f $asm_fasta contigs.fasta.sort.bam -C 2 -F 0.2 --pooled-continuous > contigs.fasta.vcf

	error_corr_assembly.pl $asm_fasta contigs.fasta.vcf > ${id}.spades.masked
	"""
}


process chewbbaca {
	publishDir "${OUTDIR}/chewbbaca", mode: 'copy', overwrite: true, pattern: '*.chewbbaca'
	cpus 7
	memory '8 GB'
	time '1h'
	// cache 'deep'
	queue='rs-fs1'

	input:
		set id, species, platform, file(asm_fasta) from maskpoly_chewbbaca

	output:
		set id, species, platform, file("${id}.chewbbaca") into chewbbaca_export
		set id, species, platform, file("${id}.missingloci") into chewbbaca_register

	script:
		cgmlst_db = params.refpath+'/species/'+species+'/cgmlst'

	"""
	flock -e /local/chewbbaca.lock \\
		  chewBBACA.py AlleleCall --fr -i ${asm_fasta} -g $cgmlst_db --ptf Staphylococcus_aureus.trn --cpu ${task.cpus} \\
		  -o chewbbaca.folder

	sed -e "s/NIPHEM/-/g" -e "s/NIPH/-/g" -e "s/LNF/-/g" -e "s/INF-*//g" -e "s/PLOT[^\t]*/-/g" -e "s/ALM/-/g" -e "s/ASM/-/g" \\
	    chewbbaca.folder/results_*/results_alleles.tsv > ${id}.chewbbaca


	tail -1 ${id}.chewbbaca | fmt -w 1 | tail -n +2 | grep '-' | wc -l > ${id}.missingloci

	"""

}


process postalignqc {
	publishDir "${OUTDIR}/postalignqc", mode: 'copy', overwrite: true
	cpus 4
	memory '8 GB'
	time '1h'

	input:
		set id, species, platform, file(bam), file(bai) from bam_postalignqc

	output:
		set id, species, platform, file("${id}.bwa.QC") into postqc_register

	script:
		cgmlst_bed = params.refpath+'/species/'+species+'/cgmlst.bed'

	"""
	postaln_qc.pl $bam $cgmlst_bed ${id} ${task.cpus} > ${id}.bwa.QC

	"""

}

process register {
	publishDir "${params.crondir}/qc", mode: 'copy', overwrite: true
	cpus 1
	memory '8 GB'
	time '1h'

	input:
		set id, species, platform, fastq_r1, fastq_r2, \
			file(quast), \
			file(postalignqc), \
			file(missingloci) \
		from fastq_register\
			 .join(quast_register,by:[0,1,2])\
			 .join(postqc_register,by:[0,1,2])\
			 .join(chewbbaca_register,by:[0,1,2])
	output:
		set id, species, platform, file("${id}.cdm")
		set id, species, platform, rundir into register_export

	script:
		parts = fastq_r1.toString().split('/')
		parts.println()
		idx =  parts.findIndexOf {it ==~ /\d{6}_.{6,8}_.{4}_.{10}/}
		rundir = parts[0..idx].join("/")
		postalignqc = params.outdir+'/'+params.subdir+'/postalignqc/'+postalignqc
		quast  = params.outdir+'/'+params.subdir+'/qc/'+quast


	"""

	read missingloci < $missingloci
	echo "--run-folder ${rundir} --sample-id ${id} --assay microbiology --qc ${postalignqc} --asmqc $quast" --micmisloc \$missingloci > ${id}.cdm

	"""
}


process to_cgviz {
	publishDir "${params.crondir}/cgviz", mode: 'copy', overwrite: true
	cpus 1
	memory '8 GB'
	time '1h'

	input:
		set id, species, platform, file(chewbbaca), \
			rundir, \
			file(quast), \
			file(mlst), \
			file(kraken), \
			file(ariba) \
		from chewbbaca_export\
			.join(register_export,by:[0,1,2])\
			.join(quast_export,by:[0,1,2])\
			.join(mlst_export,by:[0,1,2])\
			.join(kraken_export,by:[0,1,2])\
			.join(ariba_export,by:[0,1,2])

	output:
		set id, species, platform, file("${id}.cgviz")

	script:
		chewbbaca = params.outdir+'/'+params.subdir+'/chewbbaca/'+chewbbaca
		quast = params.outdir+'/'+params.subdir+'/qc/'+quast
		mlst = params.outdir+'/'+params.subdir+'/mlst/'+mlst
		kraken = params.outdir+'/'+params.subdir+'/kraken/'+kraken
		ariba = params.outdir+'/'+params.subdir+'/ariba/'+ariba

	"""
	echo "import_cgviz.pl --in $chewbbaca --overwrite --id $id --species $species --run $rundir --quast $quast --mlst $mlst --kraken $kraken --aribavir $ariba" > ${id}.cgviz
	"""
}
