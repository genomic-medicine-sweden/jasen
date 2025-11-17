#!/usr/bin/env python3

import logging
import click
import yaml

VERSION = "0.1.0"

LOG = logging.getLogger(__name__)

OUTPUT_SCHEMA_VERSION = 1

@click.group()
@click.version_option(VERSION)
@click.option("-s", "--silent", is_flag=True)
@click.option("-d", "--debug", is_flag=True)
def cli(silent, debug):
    """Jasen pipeline result processing tool."""
    if silent:
        log_level = logging.WARNING
    elif debug:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO
    # configure logging
    logging.basicConfig(
        level=log_level, format="[%(asctime)s] %(levelname)s in %(module)s: %(message)s"
    )

@cli.command()
@click.option(
    "--amrfinder",
    type=click.Path(),
    help="Amrfinderplus anti-microbial resistance results",
)
@click.option(
    "--bam",
    type=click.Path(),
    help="Read mapping to reference genome"
)
@click.option(
    "--bai",
    type=click.Path(),
    help="Read mapping index to reference genome"
)
@click.option(
    "--chewbbaca",
    type=click.Path(),
    help="cgMLST prediction results"
)
@click.option(
    "--emmtyper",
    type=click.Path(),
    help="Emmtyper m-type prediction results"
)
@click.option(
    "--gambitcore",
    type=click.Path(),
    help="Gambitcore prediction results"
)
@click.option(
    "--groups",
    multiple=True,
    required=True,
    help="Sample groups"
)
@click.option(
    "--kleborate",
    type=click.Path(),
    help="Klebsiella and Escherichia analysis pipeline"
)
@click.option(
    "--kleborate-hamronization",
    type=click.Path(),
    help="Kleborate hamronization output"
)
@click.option(
    "--kraken",
    type=click.Path(),
    help="Kraken species annotation results"
)
@click.option(
    "--lims-id",
    type=click.Path(),
    help="Lims id"
)
@click.option(
    "--mlst",
    type=click.Path(),
    help="MLST prediction results"
)
@click.option(
    "--mykrobe",
    type=click.Path(),
    help="Mykrobe results"
)
@click.option(
    "--nextflow-run-info",
    type=click.Path(),
    help="Nextflow metadata from the pipeline in json format",
)
@click.option(
    "--postalnqc",
    type=click.Path(),
    help="Postalignqc qc results"
)
@click.option(
    "--quast",
    type=click.Path(),
    help="Quast quality control metrics"
)
@click.option(
    "--ref-genome-annotation",
    type=click.Path(),
    help="Reference genome in gff format"
)
@click.option(
    "--ref-genome-sequence",
    type=click.Path(),
    help="Reference genome fasta file"
)
@click.option(
    "--resfinder",
    type=click.Path(),
    help="Resfinder resistance prediction results",
)
@click.option(
    "--sample-id",
    required=True,
    help="Sample identifier"
)
@click.option(
    "--sample-name",
    required=True,
    help="Sample name"
)
@click.option(
    "--sccmec",
    type=click.Path(),
    help="Sccmec results"
)
@click.option(
    "--serotypefinder",
    type=click.Path(),
    help="Serotypefinder serotype prediction results",
)
@click.option(
    "--shigapass",
    type=click.Path(),
    help="Shigapass results"
)
@click.option(
    "--ska-index",
    type=click.Path(),
    help="Ska index filepath"
)
@click.option(
    "--software-info",
    type=click.Path(),
    multiple=True,
    help="Software metadata from the pipeline in json format",
)
@click.option(
    "--sourmash-signature",
    type=click.Path(),
    help="Sourmash signature filepath"
)
@click.option(
    "--spatyper",
    type=click.Path(),
    help="Shigapass results"
)
@click.option(
    "--tb-grading-rules-bed",
    type=click.Path(),
    help="TB grading rules bed",
)
@click.option(
    "--tbdb-bed",
    type=click.Path(),
    help="Tbdb bed file",
)
@click.option(
    "--tbprofiler",
    type=click.Path(),
    help="Tbprofiler results"
)
@click.option(
    "--vcf",
    type=click.Path(),
    help="VCF filepath"
)
@click.option(
    "--virulencefinder",
    type=click.Path(),
    help="Virulence factor prediction results",
)
@click.option(
    "--output",
    required=True,
    type=click.Path(),
    help="Output filepath"
)
def cli(
    amrfinder,
    bam,
    bai,
    chewbbaca,
    emmtyper,
    gambitcore,
    groups,
    kleborate,
    kleborate_hamronization,
    kraken,
    lims_id,
    mlst,
    mykrobe,
    nextflow_run_info,
    postalnqc,
    quast,
    ref_genome_annotation,
    ref_genome_sequence,
    resfinder,
    sample_id,
    sample_name,
    sccmec,
    serotypefinder,
    shigapass,
    ska_index,
    software_info,
    sourmash_signature,
    spatyper,
    tb_grading_rules_bed,
    tbdb_bed,
    tbprofiler,
    vcf,
    virulencefinder,
    output
):

    prp_input = {"igv_annotations": []}

    if amrfinder:
        prp_input["amrfinder"] = amrfinder

    if bam and bai:
        prp_input["igv_annotations"].append(create_array("Read coverage", "alignment", "uri", bam, "index_uri", bai))

    if chewbbaca:
        prp_input["chewbbaca"] = chewbbaca

    if emmtyper:
        prp_input["emmtyper"] = emmtyper

    if gambitcore:
        prp_input["gambitcore"] = gambitcore

    if groups:
        prp_input["groups"] = list(groups)

    if kleborate:
        prp_input["kleborate"] = kleborate

    if kleborate_hamronization:
        prp_input["kleborate_hamronization"] = kleborate_hamronization

    if kraken:
        prp_input["kraken"] = kraken

    if lims_id:
        prp_input["lims_id"] = lims_id

    if mlst:
        prp_input["mlst"] = mlst

    if mykrobe:
        prp_input["mykrobe"] = mykrobe

    if nextflow_run_info:
        prp_input["nextflow_run_info"] = nextflow_run_info

    if postalnqc:
        prp_input["postalnqc"] = postalnqc

    if quast:
        prp_input["quast"] = quast

    if ref_genome_annotation:
        prp_input["ref_genome_annotation"] = ref_genome_annotation

    if ref_genome_sequence:
        prp_input["ref_genome_sequence"] = ref_genome_sequence

    if resfinder:
        prp_input["resfinder"] = resfinder

    if sample_id:
        prp_input["sample_id"] = sample_id

    if sample_name:
        prp_input["sample_name"] = sample_name

    if sccmec:
        prp_input["sccmec"] = sccmec

    if serotypefinder:
        prp_input["serotypefinder"] = serotypefinder

    if shigapass:
        prp_input["shigapass"] = shigapass

    if ska_index:
        prp_input["ska_index"] = ska_index

    if software_info:
        prp_input["software_info"] = list(software_info)

    if sourmash_signature:
        prp_input["sourmash_signature"] = sourmash_signature

    if spatyper:
        prp_input["spatyper"] = spatyper

    if tb_grading_rules_bed:
        prp_input["igv_annotations"].append(create_array("tbdb grading rules bed", "bed", "uri", tb_grading_rules_bed))

    if tbdb_bed:
        prp_input["igv_annotations"].append(create_array("tbdb bed", "bed", "uri", tbdb_bed))

    if tbprofiler:
        prp_input["tbprofiler"] = tbprofiler

    if vcf:
        prp_input["igv_annotations"].append(create_array("Predicted variants", "variant", "uri", vcf))

    if virulencefinder:
        prp_input["virulencefinder"] = virulencefinder

    write_yaml(output, prp_input)

def create_array(annot_name, annot_type, filepath_key, filepath_value, index_key=None, index_value=None):
    annot_array = {"name": annot_name, "type": annot_type, filepath_key: filepath_value}
    if index_key and index_value:
        annot_array[index_key] = index_value
    return annot_array

def write_yaml(output_fpath, prp_input):
    with open(output_fpath, "w") as fout:
        yaml.dump(prp_input, fout, default_flow_style=False)

if __name__ == "__main__":
    cli()
