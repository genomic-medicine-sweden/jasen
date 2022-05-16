import click
import logging
from logging.config import dictConfig
import pandas as pd
import json
from .parse import parse_resistance_pred, parse_virulence_pred, parse_qust_results, parse_cgmlst_results, parse_mlst_results, parse_species_pred
from .models.metadata import RunInformation, SoupVersion
from .models.sample import PipelineResult
from pydantic import ValidationError


dictConfig(
    {
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {
            "default": {
                "format": "[%(asctime)s] %(levelname)s in %(module)s: %(message)s",
            }
        },
        "root": {"level": "INFO"},
    }
)
LOG = logging.getLogger(__name__)

OUTPUT_SCHEMA_VERSION = 1

@click.group()
def cli():
    pass


@cli.command()
@click.option('-i', '--sample-id', required=True, help="Sample identifier")
@click.option('-u', '--run-metadata', type=click.File(), required=True, help="Analysis metadata from pipeline")
@click.option('-q', '--quast', type=click.File(), help="Quast quality control metrics")
@click.option('-p', '--process-metadata', type=click.File(), multiple=True, help="Meta data from nextflow processes")
@click.option('-k', '--kraken', type=click.File(), help="Kraken specie annotation results")
@click.option('-m', '--mlst', type=click.File(), help="MLST prediction results")
@click.option('-c', '--cgmlst', type=click.File(), help="cgMLST prediction results")
@click.option('-v', '--virulence', type=click.File(), help="Virulence factor prediction results")
@click.option('-r', '--resistance', type=click.File(), help="Resistance prediction results")
@click.argument('output', type=click.File('w'))
def create_output(sample_id, run_metadata, quast, process_metadata, kraken, mlst, cgmlst, virulence, resistance, output):
    """Combine pipeline results into a standardized json output file."""
    # base results

    run_info = RunInformation(**json.load(run_metadata))
    results = {
        'sample_id': sample_id,
        'run_metadata': {'run': run_info},
        'qc': {},
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
        results['run_metadata']['databases'] = db_info

    if quast:
        results['qc']['assembly'] = parse_qust_results(quast)
    # typing
    if mlst:
        results['mlst'] = parse_mlst_results(mlst)
    if cgmlst:
        results['cgmlst'] = parse_cgmlst_results(cgmlst)

    # resistance of different types
    if resistance:
        pred_res = json.load(resistance)
        _, amr = parse_resistance_pred(pred_res, 'amr')
        results['antimicrobial_resistance'] = amr
        results['chemical_resistance'] = parse_resistance_pred(pred_res, 'chemical')[1]
        results['environmental_resistance'] = parse_resistance_pred(pred_res, 'environmental')[1]
    # get virulence factors in sample
    if virulence:
        res = parse_virulence_pred(virulence)
        results['virulence'] = res

    if kraken:
        LOG.info('Parse kraken results')
        results['species_prediction'] = parse_species_pred(kraken)
    else:
        results['species_prediction'] = []

    try:
        output_data = PipelineResult(output_version=OUTPUT_SCHEMA_VERSION, **results)
    except ValidationError as err:
        click.secho("Input failed Validation", fg="red")
        click.secho(err)
    LOG.info(f"Storing results to: {output.name}")
    output.write(output_data.json(indent=2))
    click.secho("Finished generating pipeline output", fg="green")