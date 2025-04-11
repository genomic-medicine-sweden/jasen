process create_prp_yaml {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), val(lims_id), val(sample_name), path(nextflow_run_info), path(mykrobe), path(tbprofiler), path(bam), path(bai), path(kraken), path(postalignqc), path(quast), path(ska), path(sourmash), path(amrfinder), path(resfinder), path(resfinder_meta), path(virulencefinder), path(virulencefinder_meta), path(chewbbaca), path(emmtyper), path(mlst), path(sccmec), path(serotypefinder), path(serotypefinder_meta), path(shigapass), path(spatyper), path(vcf)
    val reference_genome
    val reference_genome_idx
    val reference_genome_gff
    val tb_grading_rules_bed
    val tbdb_bed

    output:
    tuple val(sample_id), path(output), emit: yaml

    script:
    output                          = "${sample_id}_prp.yaml"
    def access_dir                  = params.symlink_dir        ?: params.outdir
    def amrfinder_arg               = amrfinder                 ?  "--amrfinder ${params.outdir}/${params.species_dir}/amrfinderplus/${amrfinder}" : ""
    def assay_arg                   = params.assay              ?  "--assay ${params.assay}" : ""
    def bam_arg                     = bam                       ?  "--bam ${access_dir}/${params.species_dir}/${params.bam_dir}/${bam}" : ""
    def bai_arg                     = bai                       ?  "--bai ${access_dir}/${params.species_dir}/${params.bam_dir}/${bai}" : ""
    def chewbbaca_arg               = chewbbaca                 ?  "--chewbbaca ${params.outdir}/${params.species_dir}/chewbbaca/${chewbbaca}" : ""
    def emmtyper_arg                = emmtyper                  ?  "--emmtyper ${params.outdir}/${params.species_dir}/emmtyper/${emmtyper}" : ""
    def groups_arg                  = params.groups             ?  "--groups ${params.groups.join(' --groups ')}" : ""
    def kraken_arg                  = kraken                    ?  "--kraken ${params.outdir}/${params.species_dir}/kraken/${kraken}" : ""
    def lims_id_arg                 = lims_id                   ?  "--lims-id ${lims_id}" : ""
    def mlst_arg                    = mlst                      ?  "--mlst ${params.outdir}/${params.species_dir}/mlst/${mlst}" : ""
    def mykrobe_arg                 = mykrobe                   ?  "--mykrobe ${params.outdir}/${params.species_dir}/mykrobe/${mykrobe}" : ""
    def nextflow_run_info_arg       = nextflow_run_info         ?  "--nextflow-run-info ${params.outdir}/${params.species_dir}/analysis_metadata/${nextflow_run_info}" : ""
    def postalignqc_arg             = postalignqc               ?  "--postalnqc ${params.outdir}/${params.species_dir}/postalignqc/${postalignqc}" : ""
    def quast_arg                   = quast                     ?  "--quast ${params.outdir}/${params.species_dir}/quast/${quast}" : ""
    def reference_genome_arg        = reference_genome          ?  "--ref-genome-sequence ${reference_genome}" : ""
    def reference_genome_gff_arg    = reference_genome_gff      ?  "--ref-genome-annotation ${reference_genome_gff}" : ""
    def release_life_cycle_arg      = params.release_life_cycle ?  "--release-life-cycle ${params.release_life_cycle}" : ""
    def resfinder_arg               = resfinder                 ?  "--resfinder ${params.outdir}/${params.species_dir}/resfinder/${resfinder}" : ""
    def resfinder_meta_arg          = resfinder_meta            ?  "--software-info ${params.outdir}/${params.species_dir}/resfinder/${resfinder_meta}" : ""
    def sccmec_arg                  = sccmec                    ?  "--sccmec ${params.outdir}/${params.species_dir}/sccmec/${sccmec}" : ""
    def serotypefinder_arg          = serotypefinder            ?  "--serotypefinder ${params.outdir}/${params.species_dir}/serotypefinder/${serotypefinder}" : ""
    def serotypefinder_meta_arg     = serotypefinder_meta       ?  "--software-info ${params.outdir}/${params.species_dir}/serotypefinder/${serotypefinder_meta}" : ""
    def shigapass_arg               = shigapass                 ?  "--shigapass ${params.outdir}/${params.species_dir}/shigapass/${shigapass}" : ""
    def ska_dir                     = task.ext.args_ska         ?: ""
    def ska_arg                     = ska                       ?  "--ska ${ska_dir}/${ska}" : ""
    def sourmash_dir                = task.ext.args_sourmash    ?: ""
    def sourmash_arg                = sourmash                  ?  "--sourmash ${sourmash_dir}/${sourmash}" : ""
    def spatyper_arg                = spatyper                  ?  "--spatyper ${params.outdir}/${params.species_dir}/spatyper/${spatyper}" : ""
    def tb_grading_rules_bed_arg    = tb_grading_rules_bed      ?  "--tb-grading-rules-bed ${tb_grading_rules_bed}" : ""
    def tbdb_bed_arg                = tbdb_bed                  ?  "--tbdb-bed ${tbdb_bed}" : ""
    def tbprofiler_arg              = tbprofiler                ?  "--tbprofiler ${params.outdir}/${params.species_dir}/tbprofiler_mergedb/${tbprofiler}" : ""
    def vcf_arg                     = vcf                       ?  "--vcf ${access_dir}/${params.species_dir}/${params.vcf_dir}/${vcf}" : ""
    def virulencefinder_arg         = virulencefinder           ?  "--virulencefinder ${params.outdir}/${params.species_dir}/virulencefinder/${virulencefinder}" : ""
    def virulencefinder_meta_arg    = virulencefinder_meta      ?  "--software-info ${params.outdir}/${params.species_dir}/virulencefinder/${virulencefinder_meta}" : ""
    """
    create_prp_yaml.py \\
        ${amrfinder_arg} \\
        ${assay_arg} \\
        ${bam_arg} \\
        ${bai_arg} \\
        ${chewbbaca_arg} \\
        ${emmtyper_arg} \\
        ${groups_arg} \\
        ${kraken_arg} \\
        ${lims_id_arg} \\
        ${mlst_arg} \\
        ${mykrobe_arg} \\
        ${nextflow_run_info_arg} \\
        ${postalignqc_arg} \\
        ${quast_arg} \\
        ${reference_genome_arg} \\
        ${reference_genome_gff_arg} \\
        ${release_life_cycle_arg} \\
        ${resfinder_arg} \\
        ${resfinder_meta_arg} \\
        --sample-id ${sample_id} \\
        --sample-name ${sample_name} \\
        ${sccmec_arg} \\
        ${serotypefinder_arg} \\
        ${serotypefinder_meta_arg} \\
        ${shigapass_arg} \\
        ${ska_arg} \\
        ${sourmash_arg} \\
        ${spatyper_arg} \\
        ${tb_grading_rules_bed_arg} \\
        ${tbdb_bed_arg} \\
        ${tbprofiler_arg} \\
        ${vcf_arg} \\
        ${virulencefinder_arg} \\
        ${virulencefinder_meta_arg} \\
        --output ${output}
    """

    stub:
    output = "${sample_id}_prp.yaml"
    """
    touch ${output}
    """
}
