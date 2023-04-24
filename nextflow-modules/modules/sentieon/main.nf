process COPY_FASTQ {
	tag "$id"
	label "process_low"
	
	input:
		tuple	val(group),	val(id), path(read1), path(read2)

	output:
		tuple	val(group),	val(id), path("*R1*gz"), path ("*R2*gz")

	"""
	rsync -avzh ${read1} ${id}_R1.fq.gz
	rsync -avzh ${read2} ${id}_R2.fq.gz
	"""
}

process BWA_ALIGN_SHARDED {
	tag "$shard $id"
	container =   '/fs1/resources/containers/wgs_active.sif'
	label "process_high"

    input:
		val(K_size)
		val(bwa_num_shards)
        tuple	val(shard),	val(group), val(id), path(r1), path(r2) 

    output:
        tuple	val(id), path("${id}_${shard}.bwa.sort.bam"), path("${id}_${shard}.bwa.sort.bam.bai")

	when:
		params.shardbwa

	script:
	"""
	sentieon bwa mem -M \\
		-R "@RG\\tID:${id}\\tSM:${id}\\tPL:illumina" \\
		-K ${K_size} \\
		-t ${task.cpus} \\
		-p ${params.genome_file} '<sentieon fqidx extract -F ${shard}/${bwa_num_shards} -K ${K_size} ${r1} ${r2}' | sentieon util sort \\
		-r ${params.genome_file} \\
		-o ${id}_${shard}.bwa.sort.bam \\
		-t ${task.cpus} --sam2bam -i -
	"""
}

process BWA_MERGE_SHARDS {
	tag "$id"
	container =   '/fs1/resources/containers/wgs_active.sif'
	label "process_high"

	input: 
		tuple	val(id), path(shard), path(shard_bai)
	
	output:
		tuple	val(id), \
				path("${id}_merged.bam"), \
				path("${id}_merged.bam.bai")

	when:
		params.shardbwa

	script:
		bams = shard.sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() } .join(' ')

	"""
	sentieon util merge -o ${id}_merged.bam ${bams}
	"""
}

process DELETE_FASTQ {
	tag "$id"
	container = '/fs1/resources/containers/wgs_active.sif'
	label "process_low"
	
	input:
		tuple 	val(id), path(bam), path(bai)
		tuple	val(group),	val(id), path(read1), path(read2)

	output:
		tuple	val(group),	val(id)
	
	"""
	rm -r ${read1} ${read2}
	"""
}

process BAM_CRAM_ALL {
	tag "$id"
	container =   '/fs1/resources/containers/wgs_active.sif'
	label "process_high"

	input: 
		path(params.genome_file)
		tuple	val(id), path(mergedbam), path(mergedbambai)
	
	output:
		tuple	val(id), path("${id}_merged.cram"), path("${id}_merged.cram.crai") ,path("${id}_merged.cram.bai")

	script:
	"""
	sentieon driver \\
		-r ${params.genome_file} \\
		-t ${task.cpus} \\
		-i ${id}_merged.bam \\
		--algo ReadWriter \\
		--cram_write_options version=3.0  \\
		${id}_merged.cram
	
	"""
}

process LOCUS_COLLECTOR {
	tag "$shard_name $id"
	container =   '/fs1/resources/containers/wgs_active.sif'
	label "process_medium"
	maxErrors 5

	input:
		path(params.genome_file)
		tuple	val(id), path(cram), path(crai), path(bai), val(shard_name), val(shard) 

	output:
		tuple	val(id), path("${shard_name}_${id}.score"), path("${shard_name}_${id}.score.idx")
	
	"""
	sentieon driver \\
		-r ${params.genome_file} \\
		-t ${task.cpus} \\
		-i ${id}_merged.cram ${shard} \\
		--algo LocusCollector \\
		--fun score_info ${shard_name}_${id}.score
	"""
}

// This is the bottle-neck for the process completeion 
process DEDUP {
	tag "$shard_name $id"
	container =   '/fs1/resources/containers/wgs_active.sif'
	label "process_medium"

	input:
		path(params.genome_file)
		tuple	val(id), path(score), path(idx), path(cram), path(crai), path(bai), val(shard_name), val(shard)

	output:
		tuple	val(id), path("${shard_name}_${id}.deduped.cram"), path("${shard_name}_${id}.deduped.cram.crai"), path("${shard_name}_${id}.deduped.cram.bai")
		tuple	val(id), path("${shard_name}_${id}_dedup_metrics.txt")

	script:

	scores = score.sort(false) { a, b -> a.getBaseName().tokenize("_")[0] as Integer <=> b.getBaseName().tokenize("_")[0] as Integer } .join(' --score_info ')

	"""
	sentieon driver \\
		-r ${params.genome_file} \\
		-t ${task.cpus} \\
		-i  ${id}_merged.cram \\
		${shard} \\
		--algo Dedup \\
		--score_info ${scores} \\
		--rmdup  \\
		--cram_write_options version=3.0 \\
		--metrics ${shard_name}_${id}_dedup_metrics.txt \\
		${shard_name}_${id}.deduped.cram
	"""
}

process DEDUP_METRICS_MERGE {
	tag "$id"
	container = '/fs1/resources/containers/wgs_active.sif'

	input:
		tuple	val(id), path(dedup)

	output:
		tuple	val(id), path("dedup_metrics.txt") 

	"""
	sentieon driver --passthru --algo Dedup --merge dedup_metrics.txt ${dedup}
	"""
}

process SENTIEON_QC {
	tag "$id"
	container = '/fs1/resources/containers/wgs_active.sif'
	label "process_medium"
	publishDir "$params.outdir/$params.subdir/qc", 
				mode: 'copy', 
				overwrite: 'true'
	cache 'deep'
	time '2h'

	input:
		path(params.genome_file)
		tuple	val(id), path(cram), path(crai), path(bai), path(dedup)

	output:
		tuple	val(id), path("${id}.QC")

	"""
	sentieon driver \
		-r ${params.genome_file} \
		-t ${task.cpus} \
		-i ${cram} \
		--algo MeanQualityByCycle mq_metrics.txt \
		--algo QualDistribution qd_metrics.txt \
		--algo GCBias --summary gc_summary.txt gc_metrics.txt \
		--algo AlignmentStat aln_metrics.txt \
		--algo InsertSizeMetricAlgo is_metrics.txt \
		--algo WgsMetricsAlgo wgs_metrics.txt
	qc_sentieon.pl ${id} wgs > ${id}.QC
	"""
}

process BQSR {
	tag "$shard_name $id"
	container = '/fs1/resources/containers/wgs_active.sif'
	label "process_medium"
	cache 'deep'
	errorStrategy 'retry'
	maxErrors 5

	input:
		path(params.genome_file)
		path(params.KNOWN)
		tuple	val(id), path(crams), path(crai), path(bai), val(shard_name), val(shard)

	output:
		tuple	val(id), path("${shard_name}_${id}.bqsr.table")

	script:
		cram = crams.sort(false) { a, b -> a.getBaseName().tokenize("_")[0] as Integer <=> b.getBaseName().tokenize("_")[0] as Integer }
	
		Integer integerNum = shard_name as Integer
	
		if( integerNum == 1 ) {
			//def test = []
			test = cram[0..integerNum]
			cram_neigh = test.join(' -i ')
		}
		
		else if ( integerNum > 1  &&  integerNum < 31) { 
			Integer behind = integerNum - 2
			Integer forward = integerNum 
			test = cram[behind..forward]
			cram_neigh = test.join(' -i ')
		}

		else if ( integerNum == 31 ) { 
			Integer behind = integerNum - 2
			Integer forward = integerNum - 1
			test = cram[behind..forward]
			cram_neigh = test.join(' -i ')
		}

		else if ( integerNum == 32 ) { 
			Integer correction = integerNum - 1
			test = cram[correction]
			cram_neigh = test.join('')
		}
		"""
		sentieon driver \
			-r ${params.genome_file} \
			-t ${task.cpus} \
			-i $cram_neigh $shard \
			--algo QualCal \
			-k $params.KNOWN ${shard_name}_${id}.bqsr.table
		"""			
}
	
process MERGE_BQSR {
	tag "$id"
	label "process_low"
	container = '/fs1/resources/containers/wgs_active.sif'

	input:
		tuple val(id), path(tables) 

	output:
		tuple val(id), path("${id}_merged.bqsr.table") 

	"""
	sentieon driver \
		--passthru \
		--algo QualCal \
		--merge ${id}_merged.bqsr.table $tables
	"""			
}
	
process MERGE_DEDUP_CRAM {
	tag "$id"
	label "process_low"
	container = '/fs1/resources/containers/wgs_active.sif'
	publishDir "$params.outdir/$params.subdir/bam",
				mode: 'copy',
				overwrite: 'true'

	input:
		path(params.genome_file)
		tuple 	val(id), path(crams), path(crai), path(bai)

	output:
		tuple val(id), path("${id}_merged_dedup.cram"), path("${id}_merged_dedup.cram.crai"), path("${id}_merged_dedup.cram.bai")

	script:
		cram = crams.sort(false) { a, b -> a.getBaseName().tokenize("_")[0] as Integer <=> b.getBaseName().tokenize("_")[0] as Integer } .join(' -i ')
		
	"""
	sentieon util merge \
		-r ${params.genome_file} \
		-i ${cram} \
		--cram_write_options version=3.0 \
		-o ${id}_merged_dedup.cram --mergemode 10
	"""
}

process QC_TO_CDM {
	tag "$id"
	label "process_low"
	publishDir "$params.crondir/qc",
				mode: 'copy', 
				overwrite: 'true'
	
	input:
		tuple val(id), path(qc), val(diagnosis), val(r1)

	output:
		file "${id}.cdm" 

	script:
		parts = r1.toString().split("/")
		idx =  parts.findIndexOf {it ==~ /......_......_...._........../}
		rundir = parts[0..idx].join("/")

	"""
	echo "--run-folder $rundir --sample-id $id --subassay $diagnosis --assay tumwgs --qc $params.outdir/$params.subdir/qc/${id}.QC" > ${id}.cdm
	"""
}

process TNSCOPE {
	tag "$shard_name $smpl_id"
	container =   '/fs1/resources/containers/wgs_active.sif'
	label "process_medium"
	errorStrategy 'retry'
	maxErrors 5

	input:
		tuple val(shard_name), val(shard), val(group), val(smpl_id), val(type)
		tuple val(id1), path(crams1), path(crai1), path(bai1)
		tuple val(id2), path(crams2), path(crai2), path(bai2)
		tuple val(s1), path(bqsr1) 
		tuple val(s2), path(bqsr2)
		path (params.genome_file)

	output:
		tuple val("${ID_Tumor}_${ID_normal}"), path("${shard_name}.vcf.gz"), path("${shard_name}.vcf.gz.tbi")

	script:
		Integer integerNum = shard_name as Integer
		//println (shard_name)

		Tumor_index = type.findIndexOf{ it == 'tumor' || it == 'T' }
		ID_Tumor = smpl_id[Tumor_index]
		Normal_index = type.findIndexOf{ it == 'normal' || it == 'N' }
		ID_normal = smpl_id[Normal_index]

		if ( "$ID_Tumor" == id1 ) {
			cramTum = crams1.sort(false) { a, b -> a.getBaseName().tokenize("_")[0] as Integer <=> b.getBaseName().tokenize("_")[0] as Integer }
			
		}
		else if  ( "$ID_Tumor" == id2 ) {
			cramTum = crams2.sort(false) { a, b -> a.getBaseName().tokenize("_")[0] as Integer <=> b.getBaseName().tokenize("_")[0] as Integer }
		}

		if ( "$ID_normal" == id1 ) {
			cramNorm = crams1.sort(false) { a, b -> a.getBaseName().tokenize("_")[0] as Integer <=> b.getBaseName().tokenize("_")[0] as Integer }
			
		}
		else if  (  "$ID_normal" == id2 ) {
			cramNorm = crams2.sort(false) { a, b -> a.getBaseName().tokenize("_")[0] as Integer <=> b.getBaseName().tokenize("_")[0] as Integer }
		}

		if ( "$ID_Tumor" == s1 ) {
			tbqsr = bqsr1
		}
		else if ( "$ID_Tumor" == s2 ) {
			tbqsr = bqsr2
		}

		if ( "$ID_normal" == s1 ) {
			nbqsr = bqsr1
		}
		else if ( "$ID_normal" == s2 ) {
			nbqsr = bqsr2
		}

	
		if( integerNum == 1 ) {
			//def test = []
			testT = cramTum[0..integerNum]
			testN = cramNorm[0..integerNum]
			cramTum_neigh = testT.join(' -i ')
			cramNorm_neigh = testN.join(' -i ')
		}
		
		else if ( integerNum > 1  &&  integerNum < 31) { 
			Integer behind = integerNum - 2
			Integer forward = integerNum 
			testT = cramTum[behind..forward]
			testN = cramNorm[behind..forward]
			cramTum_neigh = testT.join(' -i ')
			cramNorm_neigh = testN.join(' -i ')
		}

		else if ( integerNum == 31 ) { 
			Integer behind = integerNum - 2
			Integer forward = integerNum - 1
			testT = cramTum[behind..forward]
			testN = cramNorm[behind..forward]
			cramTum_neigh = testT.join(' -i ')
			cramNorm_neigh = testN.join(' -i ')
		}

		else if ( integerNum == 32 ) { 
			Integer correction = integerNum - 1
			testT = cramTum[correction]
			testN = cramNorm[correction]
			cramTum_neigh = testT.join(' -i ')
			cramNorm_neigh = testN.join(' -i ')
		}

	    """
	    sentieon driver \
		    -t ${task.cpus} \
		    -r ${params.genome_file} \
		    -i ${cramTum_neigh} -i ${cramNorm_neigh} \
		    --cram_read_options decode_md=0 \
		    -q ${tbqsr} -q ${nbqsr} ${shard} \
		    --algo TNscope --disable_detector sv --tumor_sample ${ID_Tumor} --normal_sample ${ID_normal}  ${shard_name}.vcf.gz
	    """
}

process MERGE_VCF {
	tag "$id"
	container =   '/fs1/resources/containers/wgs_active.sif'
	label "process_medium"

	input:
		tuple	val(id), path(vcfs), path(tbi)
        
	output:
		tuple	val(group), path("${id}.tnscope.vcf.gz"), path("${id}.tnscope.vcf.gz.tbi") 

	script:
		group = "vcfs"
		vcfs_sorted = vcfs.sort(false) { a, b -> a.getBaseName().tokenize(".")[0] as Integer <=> b.getBaseName().tokenize(".") [0] as Integer } .join(' ')

	"""
	sentieon driver \\
		-t ${task.cpus} \\
		--passthru \\
		--algo TNscope \\
		--merge ${id}.tnscope.vcf.gz ${vcfs_sorted}
	"""
}

process DNASCOPE_TUM {
	tag "$group $smpl_id"
	container =   '/fs1/resources/containers/wgs_active.sif'
	label "process_high"
	errorStrategy 'retry'
	maxErrors 5

	input:
		path (params.genome_file)
		tuple val(id1), path(crams1), path(crai1), path(bai1)
		tuple val(id2), path(crams2), path(crai2), path(bai2)
		tuple val(group), val(smpl_id), val(type)

	output:
		tuple val(ID_Tumor), path("${ID_Tumor}_dnascope.vcf.gz") 
	
	script:
		Tumor_index = type.findIndexOf{ it == 'tumor' || it == 'T' }
		ID_Tumor = smpl_id[Tumor_index]
		Normal_index = type.findIndexOf{ it == 'normal' || it == 'N' }
		ID_normal = smpl_id[Normal_index]

		if ( "$ID_Tumor" == id1 ) {
			cramsTum = crams1.sort(false) { a, b -> a.getBaseName().tokenize("_")[0] as Integer <=> b.getBaseName().tokenize("_")[0] as Integer }
			cramTum = cramsTum.join(' -i ')
			}

		else if  ( "$ID_Tumor" == id2 ) {
			cramsTum = crams2.sort(false) { a, b -> a.getBaseName().tokenize("_")[0] as Integer <=> b.getBaseName().tokenize("_")[0] as Integer }
			cramTum = cramsTum.join(' -i ')

		}

		"""
		sentieon driver \
			-t ${task.cpus} \
			-r ${params.genome_file} \
			-i ${cramTum} \
			--cram_read_options decode_md=0 \
			--algo DNAscope \
			--emit_mode GVCF ${ID_Tumor}_dnascope.vcf.gz
		"""
}

process DNASCOPE_NOR {
	tag "$group $smpl_id"
	container =   '/fs1/resources/containers/wgs_active.sif'
	label "process_high"
	errorStrategy 'retry'
	maxErrors 5

	input:
		path (params.genome_file)
		tuple val(id1), path(crams1), path(crai1), path(bai1)
		tuple val(id2), path(crams2), path(crai2), path(bai2)
		tuple val(group), val(smpl_id), val(type)

	output:
		tuple val(ID_normal), path("${ID_normal}_dnascope.vcf.gz") 
	
	script:
		Tumor_index = type.findIndexOf{ it == 'tumor' || it == 'T' }
		ID_Tumor = smpl_id[Tumor_index]
		Normal_index = type.findIndexOf{ it == 'normal' || it == 'N' }
		ID_normal = smpl_id[Normal_index]

		if ( "$ID_normal" == id1 ) {
			crams = crams1.sort(false) { a, b -> a.getBaseName().tokenize("_")[0] as Integer <=> b.getBaseName().tokenize("_")[0] as Integer }
			cram = crams.join(' -i ')
			}

		else if  ( "$ID_normal" == id2 ) {
			crams = crams2.sort(false) { a, b -> a.getBaseName().tokenize("_")[0] as Integer <=> b.getBaseName().tokenize("_")[0] as Integer }
			cram = crams.join(' -i ')

		}

		"""
		sentieon driver \
			-t ${task.cpus} \
			-r ${params.genome_file} \
			-i ${cram} \
			--cram_read_options decode_md=0 \
			--algo DNAscope \
			--emit_mode GVCF ${ID_normal}_dnascope.vcf.gz
		"""
}

process FREEBAYES {
	tag "$group $id"
	label "process_low"
	errorStrategy 'retry'
	maxErrors 5

	input:
		val(mode)
		tuple val(group), val(id), path(bam), path(bai), val(bedid), path(bed),val(groupt), val(smpl_id), val(type)
		path (params.genome_file)


	output:
		tuple	val("freebayes"), val(group), path("freebayes_${bed}.vcf")
		
	script:
		if( mode == "paired" ) {
			Tumor_index = type.findIndexOf{ it == 'tumor' || it == 'T' }
			ID_Tumor = smpl_id[Tumor_index]
			tumor_index= id.findIndexOf{it == "$ID_Tumor" }
			bam_tumor = bam[tumor_index]

			Normal_index = type.findIndexOf{ it == 'normal' || it == 'N' }
			ID_normal = smpl_id[Normal_index]
			normal_index = id.findIndexOf{it == "$ID_normal" }
			bam_normal = bam[normal_index]


			"""
			freebayes -f ${params.genome_file} -t ${bed} --pooled-continuous --pooled-discrete --min-repeat-entropy 1 -F 0.03 $bam_tumor $bam_normal  > freebayes_${bed}.vcf.raw
			#vcffilter -F LowCov -f "DP > 30" -f "QA > 150" freebayes_${bed}.vcf.raw | vcffilter -F LowFrq -o -f "AB > 0.05" -f "AB = 0" | vcfglxgt > freebayes_${bed}.filt1.vcf

			filter_freebayes_somatic_wgs.pl freebayes_${bed}.vcf.raw $ID_Tumor $ID_normal | grep -v 'FAIL_' > freebayes_${bed}.vcf
			"""
		}
		else if( mode == "unpaired" ) {
			"""
			freebayes -f ${genome_file} -t ${bed} --pooled-continuous --pooled-discrete --min-repeat-entropy 1 -F 0.03 $bam > freebayes_${bed}.vcf
			"""
		}
}

process VARDICT {
	tag "$group $id"
	label "process_low"
	errorStrategy 'retry'
	maxErrors 5

	input:
		val(mode)
		tuple	val(group), val(id), path(bam), path(bai), val(bedid), path(bed), val(groupt), val(smpl_id), val(type)
		path (params.genome_file)

	output:
		tuple	val("vardict"), val(group), path("vardict_${bed}.vcf")

	script:
		if( mode == "paired" ) {
			Tumor_index = type.findIndexOf{ it == 'tumor' || it == 'T' }
			ID_Tumor = smpl_id[Tumor_index]
			tumor_index= id.findIndexOf{it == "$ID_Tumor" }
			bam_tumor = bam[tumor_index]
		
			Normal_index = type.findIndexOf{ it == 'normal' || it == 'N' }
			ID_normal = smpl_id[Normal_index]
			normal_index = id.findIndexOf{it == "$ID_normal" }
			bam_normal = bam[normal_index]
			
			"""
			export JAVA_HOME=/opt/conda/envs/CMD-TUMWGS
			vardict-java -U -th ${task.cpus} -G ${params.genome_file} -f 0.03 -N ${ID_Tumor} -b "${bam_tumor}|${bam_normal}" -c 1 -S 2 -E 3 -g 4 ${bed} | testsomatic.R | var2vcf_paired.pl -N "${ID_Tumor}|${ID_normal}" -f 0.03 > vardict_${bed}.raw.vcf
			
			filter_vardict_somatic_wgs.pl  \
				vardict_${bed}.raw.vcf ${ID_Tumor} ${ID_normal} | grep -v 'FAIL_' > vardict_${bed}.vcf

			"""
		}
		else if( mode == "unpaired" ) {
			"""
			export JAVA_HOME=/opt/conda/envs/CMD-TUMWGS
			
			vardict-java -U -G ${params.genome_file} -f 0.03 -N ${id} -b ${bam} -c 1 -S 2 -E 3 -g 4 ${bed} | teststrandbias.R | var2vcf_valid.pl -N  ${id} -E -f 0.03 > vardict_${bed}.vcf
			"""
		}

}

process CONCATENATE_VCFS {
	tag "$vc $gr"
	label "process_low"
	publishDir "$params.outdir/$params.subdir/vcf",
				mode: 'copy', 
				overwrite: true

	input:
		tuple	val(vc), val(gr), path(vcfs)
		path (params.genome_file)

	output:
		tuple	val(gr), val(vc), path("${gr}_${vc}.vcf.gz") 

	"""
	vcf-concat ${vcfs} | vcf-sort -c | gzip -c > ${vc}.concat.vcf.gz
	vt decompose ${vc}.concat.vcf.gz -o ${vc}.decomposed.vcf.gz
	vt index ${vc}.decomposed.vcf.gz
	vt sort -m chrom  ${vc}.decomposed.vcf.gz -o  ${vc}.decomposed.sorted.vcf.gz
	vt normalize ${vc}.decomposed.sorted.vcf.gz -r ${params.genome_file} | vt uniq - -o ${gr}_${vc}.vcf.gz
	"""
}

process AGGREGATE_VCFS {
	tag "$group $vc $id"
	label "process_low"
	publishDir "$params.outdir/$params.subdir/vcf",
				mode: 'copy',
				overwrite: true

	input:
		val(mode)
		tuple val(group), val(vc), path(vcfs) 
		tuple val(g), val(id), val(type) 

	output:
		tuple val(group), val("${id[tumor_idx]}"), path("${group}.agg.vcf")

	script:
	sample_order = id[0]

	if( mode == "paired" ) {
		tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }
		normal_idx = type.findIndexOf{ it == 'normal' || it == 'N' }
		sample_order = id[tumor_idx]+","+id[normal_idx]
	}

	"""
	aggregate_vcf.pl --vcf ${vcfs.sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() }.join(",")} --sample-order ${sample_order} |vcf-sort -c > ${group}.agg.unsorted.vcf
	vcf-sort -c ${group}.agg.unsorted.vcf > ${group}.agg.vcf
	"""
}

process PON_FILTER {
	tag "$group $tumor_id"
	label "process_low"
	publishDir "$params.outdir/$params.subdir/vcf", 
				mode: 'copy', 
				overwrite: true

	input:
		path (params.PON_freebayes)
		path (params.PON_vardict)
		tuple val(group), val(tumor_id), path(vcf) 

	output:
		tuple val(group), path("${group}.agg.pon.vcf")
	
	script:
		def pons = []
		pons << "freebayes="+params.PON_freebayes
		pons << "vardict="+params.PON_vardict		
		def pons_str = pons.join(",")

	"""
	filter_with_pon.pl --vcf ${vcf} --pons ${pons_str} --tumor-id $tumor_id > ${group}.agg.pon.vcf
	"""
}

process GVCF_COMBINE {
	tag "$group $id"
	label "process_medium"

	input:
		val (mode)
		tuple val(group), val(id), path(r1), path(r2)
		tuple val(id), path(vcf), path(tbi) 
	
	output:
		tuple  val(group), path("${group}.combined.vcf.gz"), path("${group}.combined.vcf.gz.tbi") 

	script:
		// Om fler än en vcf, GVCF combine annars döp om och skickade vidare
		if (mode == "family" ) {
			ggvcfs = vcf.join(' -v ')

			"""
			sentieon driver \\
				-t ${task.cpus} \\
				-r $genome_file \\
				--algo GVCFtyper \\
				-v $ggvcfs ${group}.combined.vcf.gz
			"""
		}
		// annars ensam vcf, skicka vidare
		else {
			ggvcf	= 	vcf.join('')
			gidx	= 	tbi.join('')

			"""
			mv ${ggvcf} ${group}.combined.vcf.gz
			mv ${gidx} ${group}.combined.vcf.gz.tbi
			"""
		}
}

process  SPLIT_NORMALIZE { 
	tag "$group"
	label "process_low"
	publishDir "$params.outdir/$params.subdir/vcf", 
				mode: 'copy', 
				overwrite: true

	input:
		path (params.genome_file)
		tuple	val(group), path(vcf), path(tbi)

	output:
		tuple	val(group), path("${group}.norm.uniq.DPAF.vcf") 


	"""
	vcfbreakmulti ${vcf} > ${group}.multibreak.vcf
	bcftools norm -m-both -c w -O v -f  ${params.genome_file} -o ${group}.norm.vcf ${group}.multibreak.vcf
	vcfstreamsort ${group}.norm.vcf | vcfuniq > ${group}.norm.uniq.vcf
	wgs_DPAF_filter.pl ${group}.norm.uniq.vcf > ${group}.norm.uniq.DPAF.vcf
	"""
}

process  INTERSECT {
	tag "$group"
	input:
		path (params.intersect_bed)
		tuple	 val(group), path(vcf)

	output:
		tuple val(group), path("${group}.intersected.vcf") 
		
	"""
	bedtools intersect -a ${vcf} -b ${params.intersect_bed} -u -header > ${group}.intersected.vcf
	"""
}

process ANNOTATE_VEP {
	tag "$group"
	container = '/fs1/resources/containers/ensembl-vep_latest.sif'
	label "process_high"
	publishDir  "$params.outdir/$params.subdir/vcf", 
				mode: 'copy', 
				overwrite: true

	input:
		path (params.CADD)
		path (params.VEP_FASTA)
		path (params.VEP_CACHE)
		path (params.GNOMAD)
		tuple	val(group), path(vcf)

	output:
		tuple	val(group), path("${group}.agg.pon.vep.vcf") 
		
	"""
	vep -i ${vcf} -o ${group}.agg.pon.vep.vcf \\
		--offline --merged --everything --vcf --no_stats \\
		--fork ${task.cpus} \\
		--force_overwrite \\
		--plugin CADD ${params.CADD} --plugin LoFtool \\
		--fasta ${params.VEP_FASTA} \\
		--dir_cache ${params.VEP_CACHE} \\
		--dir_plugins ${params.VEP_CACHE}/Plugins \\
		--distance 200 \\
		-cache -custom ${params.GNOMAD} \\
	"""
}

process FILTER_WITH_PANEL_SNV {
	tag "$group"
	publishDir "$params.outdir/$params.subdir/vcf",
				mode: 'copy' , 
				overwrite: true

	input:
		path	(snv_panel)
		tuple 	val(group), path(vcf) 

	output:
		 tuple 	val(group), path("${group}.agg.pon.vep.panel.vcf")

	script:
		should_hard_filter = params.SNV_HARD_FILTER ? '1' : ''

	"""
	filter_with_panel_snv.pl ${vcf} ${snv_panel} ${should_hard_filter} > ${group}.agg.pon.vep.panel.vcf

	"""
}

process MANTA {
	tag "$group $id"
	label "process_high"
	publishDir  "$params.outdir/$params.subdir/manta",
				mode:'copy', 
				overwrite: true

	input:
		val(mode)
		path(params.genome_file)
		tuple	val(group), val(smpl_id), val(type)
		tuple	val(gr), val(id), path(cram), path(crai), path(bai)

	output:
		tuple	val(group), path("${group}_manta.vcf")

	script:

	if( mode == "paired" ) {
		Tumor_index = type.findIndexOf{ it == 'tumor' || it == 'T' }
		ID_Tumor = smpl_id[Tumor_index]
		tumor_index= id.findIndexOf{it == "$ID_Tumor" }
		bam_tumor = cram[tumor_index]

		Normal_index = type.findIndexOf{ it == 'normal' || it == 'N' }
		ID_normal = smpl_id[Normal_index] 
		normal_index = id.findIndexOf{it == "$ID_normal" }
		bam_normal = cram[normal_index]

		"""
		configManta.py \\
		--tumorBam ${bam_tumor} \\
		--normalBam ${bam_normal} \\
		--reference ${params.genome_file} \\
		--runDir .

		python runWorkflow.py -m local -j ${task.cpus}

		mv ./results/variants/somaticSV.vcf.gz ${group}_manta.vcf.gz
		gunzip ${group}_manta.vcf.gz
		
		"""
		}

	else if( mode == "unpaired" ) {
		"""
		configManta.py \\
			--tumorBam ${bam} \\
			--reference ${params.genome_file} \\
			--generateEvidenceBam \\
			--region \\
			--runDir .
		
		python runWorkflow.py -m local -j ${task.cpus}
		
		mv ./results/variants/tumorSV.vcf.gz ${group}_manta.vcf.gz
		gunzip ${group}_manta.vcf.gz 
		"""
		}
}

process ANNOTATE_MANTA {
	tag "$group"
	label "process_high"
	publishDir  "$params.outdir/$params.subdir/manta" ,
				mode:'copy', 
				overwrite: true

	input:
		path(params.SNPEFF_DIR)
		tuple	val(group), path(vcf)

	output:
		tuple	val(group), path("${group}.manta.snpeff.vcf")

	"""
	snpEff -Xmx4g -configOption data.dir=${params.SNPEFF_DIR} GRCh38.86 $vcf > ${group}.manta.snpeff.vcf

	"""
}

process FILTER_WITH_PANEL_FUSIONS {
	tag "$group"
	label "process_low"
	publishDir "$params.outdir/$params.subdir/vcf",
				mode: 'copy', 
				overwrite: true

	input:
		path	(fus_panel)
		tuple	val(group), path (vcf)

	output:
		tuple	val(group), path ("${group}.manta.fusions.vcf")

	"""
	filter_with_panel_fusions.pl ${vcf} ${fus_panel} > ${group}.manta.fusions.vcf	
	"""	
}

process GATKCOV_BAF {
	tag "$group"
	label "process_high"

	input:
		path(params.GATK_GNOMAD)
		path(params.genome_file)
		tuple	val(id), val(gr), path(cram), path(crai), path(bai), val(group), val(sex), val(type)

	output:
		tuple	val(group), val(id), val(type), path ("${id}.allelicCounts.tsv")

	"""
	source activate gatk4-env
	gatk --java-options "-Xmx50g" CollectAllelicCounts \\
		-L ${params.GATK_GNOMAD} \\
		-I ${cram} \\
		-R ${params.genome_file} \\
		-O ${id}.allelicCounts.tsv
	"""	
}

process GATKCOV_COUNT_TUM {
	tag "$id $gr $sex"
	label "process_high"
	publishDir "$params.outdir/$params.subdir/cov", 
				mode: 'copy', 
				overwrite: true

	input:
		path (params.COV_INTERVAL_LIST)
		path (params.GATK_PON_FEMALE)
		path (params.GATK_PON_MALE)
		path {params.GENOMEDICT}
		path (params.genome_file)
		tuple	val(id), val(gr), path(cram), path(crai), path(bai), val(group), val(sex), val(type)

	output:
		tuple val (group), val("${id[tumor_idx]}"),  path("${id[tumor_idx]}.standardizedCR.tsv"), path("${id[tumor_idx]}.denoisedCR.tsv")
		tuple val("${id[tumor_idx]}"),  path("${id[tumor_idx]}.standardizedCR.tsv"), path("${id[tumor_idx]}.denoisedCR.tsv")

	script:
		
		PON = [F: params.GATK_PON_FEMALE, M: params.GATK_PON_MALE]
		tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T'  }

	"""
	source activate gatk4-env
	gatk CollectReadCounts \\
		-I ${cram[tumor_idx]} \\
		-L ${params.COV_INTERVAL_LIST} \\
		-R ${params.genome_file}  \\
		--interval-merging-rule OVERLAPPING_ONLY \\
		-O ${cram[tumor_idx]}.hdf5

	gatk --java-options "-Xmx50g" DenoiseReadCounts \\
		-I ${cram[tumor_idx]}.hdf5 \\
		--count-panel-of-normals ${PON[sex[tumor_idx]]} \\
		--standardized-copy-ratios ${id[tumor_idx]}.standardizedCR.tsv \\
		--denoised-copy-ratios ${id[tumor_idx]}.denoisedCR.tsv
	
	gatk PlotDenoisedCopyRatios \\
		--standardized-copy-ratios ${id[tumor_idx]}.standardizedCR.tsv \\
		--denoised-copy-ratios ${id[tumor_idx]}.denoisedCR.tsv \\
		--sequence-dictionary ${params.GENOMEDICT} \\
		--minimum-contig-length 46709983 \\
		--output . --output-prefix ${id[tumor_idx]}
    """
}

process GATKCOV_CALL_TUM {
	tag "$id $gr $sex"
	label "process_high"
	publishDir "$params.outdir/$params.subdir/cov", 
				mode: 'copy', 
				overwrite: true

	input:
		path(params.GENOMEDICT)
		tuple	val(id),  val(group), val(type), path(allelic), val(gr), path(stand), path(denoise)

	output:
		path("${id[tumor_idx]}.modeled.png")
		tuple	val("${id[tumor_idx]}"), val(group), path("${id[tumor_idx]}.called.seg")

	script:
		tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T'  }
		normal_idx = type.findIndexOf{ it == 'normal' || it == 'N'  }

	"""
	source activate gatk4-env

	gatk --java-options "-Xmx40g" ModelSegments \\
		--denoised-copy-ratios ${id[tumor_idx]}.denoisedCR.tsv \\
		--allelic-counts ${id[tumor_idx]}.allelicCounts.tsv \\
		--normal-allelic-counts ${id[normal_idx]}.allelicCounts.tsv \\
		--minimum-total-allele-count-normal 20 \\
		--output . \\
		--output-prefix ${id[tumor_idx]}

	gatk CallCopyRatioSegments \\
		--input ${id[tumor_idx]}.cr.seg \\
		--output ${id[tumor_idx]}.called.seg

	gatk PlotModeledSegments \\
		--denoised-copy-ratios ${id[tumor_idx]}.denoisedCR.tsv \\
		--allelic-counts ${id[tumor_idx]}.hets.tsv \\
		--segments ${id[tumor_idx]}.modelFinal.seg \\
		--sequence-dictionary ${params.GENOMEDICT} \\
		--minimum-contig-length 46709983 \\
		--output . \\
		--output-prefix ${id[tumor_idx]}
	"""
}

process GATKCOV_COUNT_NOR {
	tag "$id $gr $sex"
	label "process_high"
	publishDir "$params.outdir/$params.subdir/cov",	
				mode: 'copy', 
				overwrite: true

	input:
		path (params.COV_INTERVAL_LIST)
		path (params.GATK_PON_FEMALE)
		path (params.GATK_PON_MALE)
		path {params.GENOMEDICT}
		path (params.genome_file)
		tuple	val(id), val(gr), path(cram), path(crai), path(bai), val(group), val(sex), val(type)

	output:
		tuple	val(group), val("${id[normal_idx]}"),  path("${id[normal_idx]}.standardizedCR.tsv"), path("${id[normal_idx]}.denoisedCR.tsv")
		tuple	val("${id[normal_idx]}"),  path("${id[normal_idx]}.standardizedCR.tsv"), path("${id[normal_idx]}.denoisedCR.tsv")

	script:
		
		PON = [F: params.GATK_PON_FEMALE, M: params.GATK_PON_MALE]
		normal_idx = type.findIndexOf{ it == 'normal' || it == 'N'  }


	"""
	source activate gatk4-env
	
	gatk CollectReadCounts \\
		-I ${cram[normal_idx]} \\
		-L ${params.COV_INTERVAL_LIST} \\
		-R ${params.genome_file}  \\
		--interval-merging-rule OVERLAPPING_ONLY \\
		-O ${cram[normal_idx]}.hdf5

	gatk --java-options "-Xmx50g" DenoiseReadCounts \\
		-I ${cram[normal_idx]}.hdf5 \\
		--count-panel-of-normals ${PON[sex[normal_idx]]} \\
		--standardized-copy-ratios ${id[normal_idx]}.standardizedCR.tsv \\
		--denoised-copy-ratios ${id[normal_idx]}.denoisedCR.tsv
	
	gatk PlotDenoisedCopyRatios \\
		--standardized-copy-ratios ${id[normal_idx]}.standardizedCR.tsv \\
		--denoised-copy-ratios ${id[normal_idx]}.denoisedCR.tsv \\
		--sequence-dictionary ${params.GENOMEDICT} \\
		--minimum-contig-length 46709983 \\
		--output . --output-prefix ${id[normal_idx]}
    """
}

process GATKCOV_CALL_NOR {
	tag "$id $gr $sex"
	label "process_high"
	publishDir "$params.outdir/$params.subdir/cov", 
				mode: 'copy', 
				overwrite: true


	input:
		path (params.GENOMEDICT)
		tuple val(id),  val(group), val(type), path(allelic), val(gr), path(stand), path(denoise)

	output:
		path("${id[normal_idx]}.modeled.png")
		tuple  val("${id[normal_idx]}"), val(group), path("${id[normal_idx]}.called.seg")

	script:
		normal_idx = type.findIndexOf{ it == 'normal' || it == 'N'  }

	"""
	source activate gatk4-env

	gatk --java-options "-Xmx40g" ModelSegments \\
		--denoised-copy-ratios ${id[normal_idx]}.denoisedCR.tsv \\
		--output . \\
		--output-prefix ${id[normal_idx]}

	gatk CallCopyRatioSegments \\
		--input ${id[normal_idx]}.cr.seg \\
		--output ${id[normal_idx]}.called.seg

	gatk PlotModeledSegments \\
		--denoised-copy-ratios ${id[normal_idx]}.denoisedCR.tsv \\
		--segments ${id[normal_idx]}.modelFinal.seg \\
		--sequence-dictionary ${params.GENOMEDICT} \\
		--minimum-contig-length 46709983 \\
		--output . \\
		--output-prefix ${id[normal_idx]}
	"""
}

process CNVS_ANNOTATE {
	tag "$id $group"
	label "process_low"	
	publishDir "$params.outdir/$params.subdir/cnv", 
				mode: 'copy', 
				overwrite: true

	input:
		path (params.GENE_BED_PC)
		tuple  val(id), val(group), path(segments) 

	output:
		tuple	val(id), val(group), path("${id}.cnv.annotated.bed")

	"""
	overlapping_genes.pl ${segments} ${params.GENE_BED_PC} > ${id}.cnv.annotated.bed
	"""
}

process FILTER_WITH_PANEL_CNVS {
	tag "$id $group"
	label "process_low"
	publishDir "$params.outdir/$params.subdir/cnv", 
				mode: 'copy', 
				overwrite: true
	
	input:
		path	(cnv_panel)
		tuple val(id), val(group), path(bed) 

	output:
		tuple val(id), val(group), path("${id}.cnv.annotated.panel.bed") 
		
	"""
	filter_with_panel_cnv.pl ${bed} ${cnv_panel} > ${id}.cnv.annotated.panel.bed
	"""
}

process GENERATE_GENS_DATA {
	tag "$id"
	label "process_low"
	publishDir "$params.outdir/$params.subdir/plot_data", 
				mode: 'copy', 
				overwrite: 'true', 
				pattern: '*.bed.gz*'
	publishDir "$params.cron/gens", 
				mode: 'copy', 
				overwrite: 'true',
				pattern: '*.gens'			
			
	
	input:
		path(params.GENS_GNOMAD)
		tuple	val(id), path(gvcf), path(cov_stand), path(cov_denoise) 

	output:
		tuple	path("${id}.cov.bed.gz"), path("${id}.baf.bed.gz"), path("${id}.cov.bed.gz.tbi"), path("${id}.baf.bed.gz.tbi")
		path	("${id}.gens")  

	"""
	generate_gens_data.pl ${cov_stand} ${gvcf} ${id} ${params.GENS_GNOMAD}

	echo "gens load sample --sample-id $id --genome-build 38 --baf ${params.gens_accessdir}/${id}.baf.bed.gz --coverage ${params.gens_accessdir}/${id}.cov.bed.gz --overview-json ${params.gens_accessdir}/${id}.overview.json.gz" > ${id}.gens 
	"""
}

process GENERATE_GENS_DATA_NOR {
	tag "$id"
	label "process_low"
	publishDir "$params.outdir/$params.subdir/plot_data", 
				mode: 'copy', 
				overwrite: 'true', 
				pattern: '*.bed.gz*'
	publishDir "$params.cron/gens", 
				mode: 'copy', 
				overwrite: 'true',
				pattern: '*.gens'	
	input:
		path(params.GENS_GNOMAD)
		tuple	val(id), path(gvcf), path(cov_stand), path(cov_denoise) 

	output:
		tuple	path("${id}.cov.bed.gz"), path("${id}.baf.bed.gz"), path("${id}.cov.bed.gz.tbi"), path("${id}.baf.bed.gz.tbi")
		path 	("${id}.gens")

	"""
	generate_gens_data.pl ${cov_stand} ${gvcf} ${id} ${params.GENS_GNOMAD}
	
	echo "gens load sample --sample-id $id --genome-build 38 --baf ${params.gens_accessdir}/${id}.baf.bed.gz --coverage ${params.gens_accessdir}/${id}.cov.bed.gz --overview-json ${params.gens_accessdir}/${id}.overview.json.gz" > ${id}.gens 
	"""
}

process COYOTE {
	tag "$id"
	label "process_low"	
	publishDir "${params.crondir}/coyote", 
				mode: 'copy', 
				overwrite: true

	input:
		tuple 	val(group), path(vcf), val(id), path(cnv), path(fusions)
		tuple	val(g), val(type), val(lims_id), val(pool_id) 
		path (cnvplot)

	output:
		path("${group}.coyote_wgs")

	script:
		tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }

	"""
	echo "import_myeloid_to_coyote_vep_gms_dev_WGS.pl \\
		--id ${group} --group tumwgs \\
		--vcf /access/tumwgs/vcf/${vcf} \\
		--cnv /access/tumwgs/cnv/${cnv} \\
		--transloc /access/tumwgs/vcf/${fusions} \\
		--cnvprofile /access/tumwgs/cov/${cnvplot} \\
		--clarity-sample-id ${lims_id[tumor_idx]} \\
		--build 38 \\
        --gens ${group} \\
		--clarity-pool-id ${pool_id[tumor_idx]}" > ${group}.coyote_wgs
	"""
}