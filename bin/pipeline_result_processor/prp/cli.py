import click
import logging
from logging.config import dictConfig
import pandas as pd
import json
from .parse import parse_resistance_pred, parse_virulence_pred, parse_qust_results, parse_cgmlst_results, parse_mlst_results

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


SPP_MIN_READ_FRAC = 0.001


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
    results = {
        'sample_id': sample_id,
        'meta': {'run': json.load(run_metadata)},
        'qc': {},
    }
    import pdb; pdb.set_trace()
    if process_metadata:
        results['meta']['databases'] = [json.load(jf) for jf in process_metadata]

    if quast:
        results['qc']['assembly'] = dc_asdict(parse_qust_results(quast))
    # typing
    if mlst:
        results['mlst'] = dc_asdict(parse_mlst_results(mlst))
    if cgmlst:
        results['cgmlst'] = dc_asdict(parse_cgmlst_results(cgmlst))

    # resistance of different types
    if resistance:
        pred_res = json.load(resistance)
        _, amr = parse_resistance_pred(pred_res, 'amr')
        results['antimicrobial_resistance'] = amr
        results['chemical_resistance'] = parse_resistance_pred(pred_res, 'chemical')[1]
        results['environmental_resistance'] = parse_resistance_pred(pred_res, 'environmental')[1]
    # get virulence factors in sample
    if virulence:
        results['virulence'] = parse_virulence_pred(virulence)

    if kraken:
        LOG.info('Parse kraken results')
        specie_pred = pd.read_csv(kraken, sep='\t').sort_values('fraction_total_reads', ascending=False)
        specie_pred = specie_pred[specie_pred['fraction_total_reads'] > SPP_MIN_READ_FRAC]
        # limit the number of predicted species
        results['species_prediction'] = specie_pred.to_dict(orient='records')

    LOG.info(f"Storing results to: {output.name}")
    json.dump(results, output, indent=2)