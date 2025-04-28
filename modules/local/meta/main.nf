process save_analysis_metadata {
    tag "${sample_id}"
    scratch params.scratch

    input:
    tuple val(sample_id), path(reads), val(sequencing_run), val(lims_id), val(sample_name)
    val assay
    val platform
    val release_life_cycle

    output:
    tuple val(sample_id), path(output), emit: json

    script:
    def assay = assay ?: ""
    def config_files = workflow.configFiles.collect { "\"${it.toString().trim()}\"" }.unique().join(", ")
    def lims_id = lims_id ?: ""
    def profiles = workflow.profile.split(",").collect { "\"${it.trim()}\"" }.join(", ")
    def release_life_cycle = release_life_cycle ?: ""
    def sequencing_type = reads.size() == 2 ? "PE" : "SE"
    def sequencing_run = sequencing_run ?: ""
    output = "${sample_id}_analysis_meta.json"
    """
    #!/usr/bin/env python
    import json

    res = {
        "workflow_name": "${workflow.runName}",
        "sample_name": "${sample_name}",
        "lims_id": "${lims_id}",
        "assay": "${assay}",
        "release_life_cycle": "${release_life_cycle}",
        "sequencing_run": "${sequencing_run}",
        "sequencing_platform": "${platform}",
        "sequencing_type": "${sequencing_type}",
        "date": "${workflow.start}",
        "pipeline": "${workflow.scriptName}",
        "version": "${workflow.manifest.version}",
        "commit": "${workflow.commitId}",
        "configuration_files": list([${config_files}]),
        "analysis_profile": list([${profiles}]),
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
