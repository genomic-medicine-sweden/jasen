process save_analysis_metadata {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads), val(sequencing_run), val(lims_id), val(sampleName)
    val platform

    output:
    tuple val(sample_id), path(output), emit: json

    script:
    def sequencing_type = reads.size() == 2 ? "PE" : "SE"
    sequencing_run = sequencing_run ? "${sequencing_run}" : ""
    lims_id = lims_id ? "${lims_id}" : ""
    output = "${sample_id}_analysis_meta.json"
    """
    #!/usr/bin/env python
    import json

    res = {
        "workflow_name": "${workflow.runName}",
        "sample_name": "${sampleName}",
        "lims_id": "${lims_id}",
        "sequencing_run": "${sequencing_run}",
        "sequencing_platform": "${platform}",
        "sequencing_type": "${sequencing_type}",
        "date": "${workflow.start}",
        "pipeline": "${workflow.scriptName}",
        "version": "${workflow.manifest.version}",
        "commit": "${workflow.commitId}",
        "configuration_files": "${workflow.configFiles}"[1:-1].split(','),
        "analysis_profile": "${workflow.profile}",
        "command": "${workflow.commandLine}",
    }
    with open("${output}", 'w') as out:
        json.dump(res, out, indent=2)
    """

    stub:
    output = "${sample_id}_analysis_meta.json"
    """
    touch ${output}
    """
}
