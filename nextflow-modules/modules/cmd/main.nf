process export_to_cdm {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(quast), path(postQc), path(cgmlst), val(rundir)
    val species

  output:
    tuple val(sampleName), path(output), emit: cdm

  script:
    output = "${sampleName}.cdm"
    rundir = rundir ? rundir : "testdir"
    """
    cgmlstMissingLoci=\$(grep -n -o -E 'LNF|PLOT3|PLOT5|NIPH|NIPHEM|ALM|ASM' ${cgmlst} | cut -d : -f 1 | uniq -c | cut -d" " -f5)

    echo --run-folder ${rundir} \\
         --sample-id ${sampleName} \\
         --assay ${species} \\
         --qc ${postQc} \\
         --asmqc ${quast} \\
         --micmisloc \$cgmlstMissingLoci > ${output}
    """

  stub:
    output = "${sampleName}.cdm"
    """
    touch $output
    """
}
