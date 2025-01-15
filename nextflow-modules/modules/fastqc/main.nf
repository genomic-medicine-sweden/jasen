process fastqc {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(reads)

  output:
    tuple val(sampleID), path(summary_output) , emit: summary
    tuple val(sampleID), path(output)         , emit: output
    path "*versions.yml"                      , emit: versions

  script:
    def reads_arg = reads.size() == 2 ? "${reads[0]} ${reads[1]}" : "${reads[0]}"
    output_dir = "fastqc_outdir"
    summary_output = "${sampleID}_fastqc_summary.txt"
    output = "${sampleID}_fastqc_data.txt"
    """
    fastqc ${args} --extract -o ${output_dir} ${reads_arg}
    mv ${output_dir}/*/summary.txt ${summary_output}
    mv ${output_dir}/*/fastqc_data.txt ${output}

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     fastqc:
      version: \$(echo \$(fastqc --version 2>&1) | sed -r 's/^.*FastQC v//')
      container: ${task.container}
    END_VERSIONS
    """

 stub:
    """
    mkdir -p output_dir/${sampleID}
    touch output_dir/${sampleID}/summary.txt
    touch output_dir/${sampleID}/fastqc_data.txt

    cat <<-END_VERSIONS > ${sampleID}_${task.process}_versions.yml
    ${task.process}:
     fastqc:
      version: \$(echo \$(fastqc --version 2>&1) | sed -r 's/^.*FastQC v//')
      container: ${task.container}
    END_VERSIONS
    """
}
