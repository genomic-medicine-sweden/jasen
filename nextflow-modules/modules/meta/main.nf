process save_analysis_metadata {
  tag "${sampleID}"
  scratch params.scratch

  input:
    tuple val(sampleID), path(reads), val(platform), val(sequencingRun), val(limsID), val(sampleName)

  output:
    tuple val(sampleID), path(output), emit: meta

  script:
    def seqType = reads.size() == 2 ? "PE" : "SE"
    sequencingRun = sequencingRun ? "${sequencingRun}" : ""
    limsID = limsID ? "${limsID}" : ""
    output = "${sampleID}_analysis_meta.json"
    """
    #!/usr/bin/env python
    import json

    res = {
        "workflow_name": "$workflow.runName",
        "sample_name": "${sampleName}",
        "lims_id": "${limsID}",
        "sequencing_run": "${sequencingRun}",
        "sequencing_platform": "${platform}",
        "sequencing_type": "${seqType}",
        "date": "$workflow.start",
        "pipeline": "$workflow.scriptName",
        "version": "$workflow.manifest.version",
        "commit": "$workflow.commitId",
        "configuration_files": "$workflow.configFiles"[1:-1].split(','),
        "analysis_profile": "$workflow.profile",
        "command": "$workflow.commandLine",
    }
    with open("$output", 'w') as out:
        json.dump(res, out, indent=2)
    """

  stub:
    output = "${sampleID}_analysis_meta.json"
    """
    touch $output
    """
}
