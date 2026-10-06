process sccmec {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path(output), emit: tsv 
    tuple val("${task.process}"), val('sccmec'), eval("echo \$(sccmec --version 2>&1) | sed -n 's/.*sccmec_targets, version //p' | sed 's/ .*//'"), topic: versions, emit: versions_sccmec

    when:
    task.ext.when

    script:
    def args = task.ext.args ?: ''
    output = "${sample_id}_sccmec.tsv"
    outputDir = "sccmec_outdir"
    """
    sccmec --input ${assembly} --prefix ${sample_id}_sccmec -o ${outputDir} ${args}
    cp ${outputDir}/${sample_id}_sccmec.tsv ${output}
    """

    stub:
    output = "${sample_id}_sccmec.tsv"
    """
    touch ${output}
    """
}
