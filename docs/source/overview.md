# Pipeline overview

Jasen is a pipeline for resistance, virulence and epityping of whole genome sequenced infectious bacteria. The pipelines is intended to be used in clinical routine and therefore contains a currated, and opinionated, set of softwares and methods for assessing the samples quality.

The pipeline is divided into workflows dedicated for analysing one or a set of related species, for example can the *Escherichia coli* workflow also be used for analysing Shigella. See the documentation on the individual workflow for information on how it can be applied.

```{figure} _static/flowchart.png
:alt: jasen flowchart

Flowchart depiciting the main workflows, its output, and processes the workflow are constituted of. The colour of the line denotes the workflow and the area with yellow background indicate quality control processes.
```

## Status on workflow development

Jasen is being continiously developed with the aim to add support for more species. Therefore can workflows be at different stages of development and evolve at different pace. A workflow can be either of three stages (*development*, *draft*, and *stable*) of development, these reflect the stability and the likelyhood of radical canges being made to the workflow.

Workflows under **development** are being actively worked at. While all workflows should be functional in the *masters* branch and for every release, there are no guarantee that workflows under development produce accurate results and major changes can be introduce without notice. Use a your own descression.

**Draft** workflows are considered feature complete and in the process of being validated. Software versions can still change without notice. Use with descression.

**Stable** workflows are feature complete and have been validated for introduction into clinical routine by one or more partners. Updates to these workflows, and versions of softwares used by these workflows, are introduced periodically after a discussion in the Jasen community.

## Workflows

| Species                                                                   | Development status |
|---------------------------------------------------------------------------|--------------------|
| [*Staphylococcus arueus*](./workflows/staphylococcus_aureus.md)           | Draft              |
| [*Escherichia coli*](./workflows/escherichia_coli.md)                     | Draft              |
| [*Mycobacterium tuberculosis*](./workflows/mycobacterium_tuberculosis.md) | Draft              |

## Input

Jasen is designed for analysis of sequenced whole genomes from bacterial isolates. The pipeline takes short-read sequenced reads, in FastQ format, from Illumina and Ion Torrent as input.

The paths to the reads are defined in a CSV file with sample id, and optionally sequencing run and lims id. See [the usage](./useage.md) page for information on how to start a analysis.

## Output

Jasen will publish the analysis result to the path specified by the `outdir` variable in the config file. The pipeline will combine the output of the different softwares into standardised result file in JSON format for easier downstream processing of the result. The combined output have the same format regardless for workflow and can be uploaded to the result visualisation tool Bonsai for easy analysis. 

The pipeline will also publish the output files of every tool beign run seperatly in the `outdir` folder.