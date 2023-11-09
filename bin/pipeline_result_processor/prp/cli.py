import json
import logging

import click
from typing import List
from pydantic import ValidationError

from .models.metadata import RunInformation, SoupVersion
from .models.phenotype import ElementType
from .models.typing import TypingMethod
from .models.qc import QcMethodIndex
from .models.sample import MethodIndex, PipelineResult
from .parse import (
    parse_cgmlst_results,
    parse_mlst_results,
    parse_quast_results,
    parse_postalignqc_results,
    parse_resfinder_amr_pred,
    parse_amrfinder_amr_pred,
    parse_kraken_result,
    parse_virulencefinder_vir_pred,
    parse_amrfinder_vir_pred,
    parse_mykrobe_amr_pred,
    parse_mykrobe_lineage_results,
    parse_tbprofiler_amr_pred,
    parse_tbprofiler_lineage_results,
)

logging.basicConfig(
    level=logging.INFO, format="[%(asctime)s] %(levelname)s in %(module)s: %(message)s"
)
LOG = logging.getLogger(__name__)

OUTPUT_SCHEMA_VERSION = 1


@click.group()
def cli():
    pass


@cli.command()
@click.option("-i", "--sample-id", required=True, help="Sample identifier")
@click.option(
    "-u",
    "--run-metadata",
    type=click.File(),
    required=True,
    help="Analysis metadata from the pipeline in json format",
)
@click.option("-q", "--quast", type=click.File(), help="Quast quality control metrics")
@click.option(
    "-d",
    "--process-metadata",
    type=click.File(),
    multiple=True,
    help="Nextflow processes metadata from the pipeline in json format",
)
@click.option(
    "-k", "--kraken", type=click.File(), help="Kraken species annotation results"
)
@click.option(
    "-a", "--amr", type=str, help="amrfinderplus anti-microbial resistance results"
)
@click.option("-m", "--mlst", type=click.File(), help="MLST prediction results")
@click.option("-c", "--cgmlst", type=click.File(), help="cgMLST prediction results")
@click.option("-v", "--virulence", type=click.File(), help="Virulence factor prediction results")
@click.option("-r", "--resistance", type=click.File(), help="resfinder resistance prediction results")
@click.option("-p", "--quality", type=click.File(), help="postalignqc qc results")
@click.option("-k", "--mykrobe", type=click.File(), help="mykrobe results")
@click.option("-t", "--tbprofiler", type=click.File(), help="tbprofiler results")
@click.option("--correct_alleles", is_flag=True, help="Correct alleles")
@click.argument("output", type=click.File("w"))
def create_output(
    sample_id,
    run_metadata,
    quast,
    process_metadata,
    kraken,
    mlst,
    cgmlst,
    virulence,
    amr,
    resistance,
    quality,
    mykrobe,
    tbprofiler,
    correct_alleles,
    output,
):
    """Combine pipeline results into a standardized json output file."""
    # base results
    LOG.info("Start generating pipeline result json")

    run_info = RunInformation(**json.load(run_metadata))
    results = {
        "sample_id": sample_id,
        "run_metadata": {"run": run_info},
        "qc": [],
        "typing_result": [],
        "element_type_result": [],
    }
    if process_metadata:
        db_info: List[SoupVersion] = []
        for soup in process_metadata:
            dbs = json.load(soup)
            if isinstance(dbs, (list, tuple)):
                for db in dbs:
                    db_info.append(SoupVersion(**db))
            else:
                db_info.append(SoupVersion(**dbs))
        results["run_metadata"]["databases"] = db_info

    # qc
    if quast:
        LOG.info("Parse quast results")
        res: QcMethodIndex = parse_quast_results(quast)
        results["qc"].append(res)
    if quality:
        LOG.info("Parse quality results")
        res: QcMethodIndex = parse_postalignqc_results(quality)
        results["qc"].append(res)

    # typing
    if mlst:
        LOG.info("Parse mlst results")
        res: MethodIndex = parse_mlst_results(mlst)
        results["typing_result"].append(res)
    if cgmlst:
        LOG.info("Parse cgmlst results")
        res: MethodIndex = parse_cgmlst_results(cgmlst, correct_alleles=correct_alleles)
        results["typing_result"].append(res)

    # resfinder of different types
    if resistance:
        LOG.info("Parse resistance results")
        pred_res = json.load(resistance)
        methods = [ElementType.AMR, ElementType.BIOCIDE, ElementType.HEAT]
        for method in methods:
            res: MethodIndex = parse_resfinder_amr_pred(pred_res, method)
            # exclude empty results from output
            if len(res.result.genes) > 0 and len(res.result.mutations) > 0:
                results["element_type_result"].append(res)

    # amrfinder
    if amr:
        LOG.info("Parse amr results")
        methods = [
            ElementType.AMR,
            ElementType.BIOCIDE,
            ElementType.METAL,
            ElementType.HEAT,
        ]
        for method in methods:
            res: MethodIndex = parse_amrfinder_amr_pred(amr, method)
            results["element_type_result"].append(res)
        vir: MethodIndex = parse_amrfinder_vir_pred(amr)
        results["element_type_result"].append(vir)

    # get virulence factors in sample
    if virulence:
        LOG.info("Parse virulence results")
        vir: MethodIndex = parse_virulencefinder_vir_pred(virulence)
        results["element_type_result"].append(vir)

    # species id
    if kraken:
        LOG.info("Parse kraken results")
        results["species_prediction"] = parse_kraken_result(kraken)
    else:
        results["species_prediction"] = []

    # mycobacterium tuberculosis
    # mykrobe
    if mykrobe:
        LOG.info("Parse mykrobe results")
        pred_res = json.load(mykrobe)
        sample_id = list(pred_res.keys())[0]
        db_info: List[SoupVersion] = []
        db_info = [
            SoupVersion(
                **{
                    "name": "mykrobe-predictor",
                    "version": pred_res[sample_id]["version"]["mykrobe-predictor"],
                    "type": "database",
                }
            )
        ]
        results["run_metadata"]["databases"] = db_info
        amr_res: MethodIndex = parse_mykrobe_amr_pred(pred_res[sample_id], ElementType.AMR)
        results["element_type_result"].append(amr_res)
        lin_res: MethodIndex = parse_mykrobe_lineage_results(pred_res[sample_id], TypingMethod.LINEAGE)
        results["typing_result"].append(lin_res)

    # tbprofiler
    if tbprofiler:
        LOG.info("Parse tbprofiler results")
        pred_res = json.load(tbprofiler)
        db_info: List[SoupVersion] = []
        db_info = [
            SoupVersion(
                **{
                    "name": pred_res["db_version"]["name"],
                    "version": pred_res["db_version"]["commit"],
                    "type": "database",
                }
            )
        ]
        results["run_metadata"]["databases"].extend(db_info)
        lin_res: MethodIndex = parse_tbprofiler_lineage_results(pred_res, TypingMethod.LINEAGE)
        results["typing_result"].append(lin_res)
        amr_res: MethodIndex = parse_tbprofiler_amr_pred(pred_res, ElementType.AMR)
        results["element_type_result"].append(amr_res)

    try:
        output_data = PipelineResult(schema_version=OUTPUT_SCHEMA_VERSION, **results)
    except ValidationError as err:
        click.secho("Input failed Validation", fg="red")
        click.secho(err)
        raise click.Abort
    LOG.info(f"Storing results to: {output.name}")
    output.write(output_data.json(indent=2))
    click.secho("Finished generating pipeline output", fg="green")


@cli.command()
@click.argument("output", type=click.File("w"), default="-")
def print_schema(output):
    """Print Pipeline result output format schema."""
    click.secho(PipelineResult.schema_json(indent=2))


@cli.command()
@click.argument("output", type=click.File("r"))
def validate(output):
    """Validate output format of result json file."""
    js = json.load(output)
    try:
        PipelineResult(**js)
    except ValidationError as err:
        click.secho("Invalid file format X", fg="red")
        click.secho(err)
    else:
        click.secho(f'The file "{output.name}" is valid', fg="green")
