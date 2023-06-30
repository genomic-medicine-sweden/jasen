process export_to_cdm {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(cgmlstMissingLoci), path(quast), path(postQc)

  output:
    path(output)

  when:
    task.ext.when

  script:
    output = "${sampleName}.cdm"
    rundir = 'fool'

    """
    echo --run-folder ${rundir} \\
         --sample-id ${sampleName} \\
         --assay microbiology \\
         --qc ${postQc} \\
         --asmqc ${quast} \\
         --micmisloc ${cgmlstMissingLoci} > ${output}
    """

  stub:
    output = "${sampleName}.cdm"
    """
    touch $output
    """
}
