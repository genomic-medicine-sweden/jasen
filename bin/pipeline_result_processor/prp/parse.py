
#! /usr/bin/env python
"""Combine output files from pipeline into a standardized json output."""

import click
import json
import logging
from logging.config import dictConfig
import csv
import pandas as pd
from typing import Dict, Any, List, Optional
from collections import defaultdict
from dataclasses import dataclass
from dataclasses import asdict as dc_asdict
from .models.sample import AssemblyQc

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


@dataclass
class DatabaseReference:
    ref_database: str
    ref_id: str

@dataclass
class GeneBase:
    """Container for gene information"""
    name: str
    accession: str
    # prediction info
    depth: Optional[float]
    identity: float
    coverage: float
    ref_start_pos: int
    ref_end_pos: int
    ref_gene_length: int
    alignment_length: int


@dataclass
class ResistanceGene(GeneBase, DatabaseReference):
    """Container for resistance gene information"""
    phenotypes: List[str]


@dataclass
class VirulenceGene(GeneBase, DatabaseReference):
    """Container for virulence gene information"""
    virulence_category: str


@dataclass
class VariantBase:
    """Container for mutation information"""
    variant_type: str  # type of mutation insertion/deletion/substitution
    genes: List[str]
    position: int
    ref_codon: str
    alt_codon: str
    # prediction info
    depth: float


@dataclass
class ResistanceVariant(VariantBase, DatabaseReference):
    """Container for resistance variant information"""
    phenotypes: List[str]


@dataclass
class MlstResult:
    """Container for MLST results."""
    sequence_type: Optional[int]
    scheme: str  # sceme name
    alleles: Dict[str, Optional[int]]


@dataclass
class CgMlstResult:
    """Container for MLST results."""
    n_novel: int
    n_missing: int
    alleles: Dict[str, Optional[int]]


def _normalize_names(name: str) -> str:
    """Normalize names to be easier to computationally process.

    - space to underline
    - ensure all lowercase
    """
    return name.replace(' ', '_').lower()


def parse_qust_results(file: str) -> AssemblyQc:
    """Parse quast file and extract relevant metrics.

    Args:
        sep (str): seperator

    Returns:
        AssemblyQc: list of key-value pairs
    """
    LOG.info(f"Parsing tsv file: {file}")
    creader = csv.reader(file, delimiter='\t')
    header = next(creader)
    raw = [dict(zip(header, row)) for row in creader]
    return AssemblyQc(
            total_length=int(raw[0]["Total length"]),
            reference_length=raw[0]["Reference length"],
            largest_contig=raw[0]["Largest contig"],
            n_contigs=raw[0]["# contigs"],
            n50=raw[0]["N50"],
            assembly_gc=raw[0]["GC (%)"],
            reference_gc=raw[0]["Reference GC (%)"],
            duplication_ratio=raw[0]["Duplication ratio"]
    )


def parse_mlst_results(path: str) -> MlstResult:
    """Parse mlst results from mlst to json object."""
    LOG.info("Parsing mlst results")
    result = json.load(path)[0]
    result_obj = MlstResult(
        scheme=result['scheme'],
        sequence_type=None if result['sequence_type'] == '-' else result['sequence_type'],
        alleles={
            gene: None if allele == '-' else allele
            for gene, allele
            in result['alleles'].items()
        },
    )
    return result_obj


def parse_cgmlst_results(file: str, include_novel_alleles: bool = True) -> Dict[str, str]:
    """Parse chewbbaca cgmlst prediction results to json results.

    chewbbaca reports errors in allele profile, https://github.com/B-UMMI/chewBBACA
    -------------------
    INF-<allele name>, inferred new allele
    LNF, loci not found
    PLOT, loci contig tips
    NIPH, non-informative paralogous hits
    NIPHEM,
    ALM, alleles larger than locus length
    ASM, alleles smaller than locus length
    """
    ERRORS = ['LNF', 'PLOT5', 'NIPH', 'NIPHEM', 'ALM', 'ASM']

    def replace_errors(allele):
        """Replace errors and novel alleles with nulls if they are not to be inlcuded."""
        if any([allele in ERRORS, allele.startswith("INF") and not include_novel_alleles]):
            return None
        elif allele.startswith("INF") and include_novel_alleles:
            return int(allele.split('-')[1])
        return int(allele)

    msg = "Parsing cgmslt results, "
    LOG.info(msg + "not" if not include_novel_alleles else "" + "including novel alleles")
    creader = csv.reader(file, delimiter='\t')
    _, *allele_names = (colname.rstrip('.fasta') for colname in next(creader))
    # parse alleles
    _, *alleles = next(creader)
    corrected_alleles = (replace_errors(a) for a in alleles)
    results = CgMlstResult(
            n_novel=sum(1 for a in alleles if a.startswith('INF')),
            n_missing=sum(1 for a in alleles if a in ERRORS),
            alleles=dict(zip(allele_names, corrected_alleles)),
    )
    return results


def  _get_resfinder_amr_sr_profie(resfinder_result, limit_to_phenotypes=None):
    """Get resfinder susceptibility/resistance profile."""
    susceptible = set()
    resistant = set()
    for phenotype in resfinder_result['phenotypes'].values():
        # skip phenotype if its not part of the desired category
        if limit_to_phenotypes is not None and phenotype['key'] not in limit_to_phenotypes:
            continue

        if phenotype['resistant']:
            resistant.add(phenotype['resistance'])
        else:
            susceptible.add(phenotype['resistance'])
    return {'susceptible': list(susceptible), 'resistant': list(resistant)}


def _parse_resfinder_amr_genes(resfinder_result, limit_to_phenotypes=None):
    """Get resistance genes from resfinder result."""
    results = []
    for info in resfinder_result['genes'].values():
        # Get only acquired resistance genes
        if not info['ref_database'].startswith('Res'):
            continue

        # Get only genes of desired phenotype
        if limit_to_phenotypes is not None:
            intersect = set(info['phenotypes']) & set(limit_to_phenotypes)
            if len(intersect) == 0:
                continue

        # store results
        gene = ResistanceGene(
            name=info['name'],
            accession=info['ref_acc'],
            depth=info['depth'],
            identity=info['identity'],
            coverage=info['coverage'],
            ref_start_pos=info['ref_start_pos'],
            ref_end_pos=info['ref_end_pos'],
            ref_gene_length=info['ref_gene_lenght'],
            alignment_length=info['alignment_length'],
            phenotypes=info['phenotypes'],
            ref_database=info['ref_database'],
            ref_id=info['ref_id']
        )
        results.append(dc_asdict(gene))
    return results


def _parse_resfinder_amr_variants(resfinder_result, limit_to_phenotypes=None):
    """Get resistance genes from resfinder result."""
    results = []
    for info in resfinder_result['seq_variations'].values():
        # Get only variants from desired phenotypes
        if limit_to_phenotypes is not None:
            intersect = set(info['phenotypes']) & set(limit_to_phenotypes)
            if len(intersect) == 0:
                continue
        # get gene depth
        info['depth'] = resfinder_result['genes'][info['genes'][0]]['depth']
        # translate variation type bools into classifier
        if info['substitution']:
            var_type = 'substitution'
        elif info['insertion']:
            var_type = 'insertion'
        elif info['deletion']:
            var_type = 'deletion'
        else:
            raise ValueError('Output has no known mutation type')
        variant = ResistanceVariant(
            variant_type=var_type,
            genes=info['genes'],
            phenotypes=info['phenotypes'],
            position=info['ref_start_pos'],
            ref_codon=info['ref_codon'],
            alt_codon=info['var_codon'],
            depth=info['depth'],
            ref_database=info['ref_database'],
            ref_id=info['ref_id'],
        )
        results.append(dc_asdict(variant))
    return results


def parse_resistance_pred(prediction: Dict[str, Any], resistance_category) -> Dict[str, str]:
    """Parse resfinder resistance prediction results."""
    # resfinder missclassifies resistance the param amr_category by setting all to amr
    LOG.info("Parsing resistance prediction")
    meta = [{
        'name': metrics['database_name'],
        'version': metrics['database_version'],
        'type': metrics['type'],
        }
        for metrics
        in prediction['databases'].values()
    ]
    # parse resistance based on the category
    categories = {
        'chemical': ['formaldehyde',  'benzylkonium chloride', 'ethidium bromide', 'chlorhexidine', 'cetylpyridinium chloride', 'hydrogen peroxide'],
        'environmental': ['temperature']
    }
    categories['amr'] = list({k for k in prediction['phenotypes'].keys()} - set(categories['chemical'] + categories['environmental']))

    # parse resistance
    resistance = {
        'prediction': _get_resfinder_amr_sr_profie(prediction, categories[resistance_category]),
        'genes': _parse_resfinder_amr_genes(prediction, categories[resistance_category]),
        'mutations': _parse_resfinder_amr_variants(prediction, categories[resistance_category]),
    }
    return meta, resistance


def parse_virulence_pred(file: str) -> Dict[str, str]:
    """Parse virulence prediction results."""
    LOG.info("Parsing virulence prediction")
    pred = json.load(file)
    results = {}
    if 'virulencefinder' in pred:
        # parse vir[bb[jjjjulence finder results
        species = [k for k in pred['virulencefinder']['results'].keys()]
        for key, genes in pred['virulencefinder']['results'][species[0]].items():
            virulence_category = key.split('_')[1]
            vir_genes = []
            for gn in genes.values():
                start_pos, end_pos = map(int, gn['position_in_ref'].split('..'))
                gene = VirulenceGene(
                        name=gn['virulence_gene'],
                        virulence_category=virulence_category,
                        accession=gn['accession'],
                        depth=None,
                        identity=gn['identity'],
                        coverage=gn['coverage'],
                        ref_start_pos=start_pos,
                        ref_end_pos=end_pos,
                        ref_gene_length=gn['template_length'],
                        alignment_length=gn['HSP_length'],
                        ref_database='virulenceFinder',
                        ref_id=gn['hit_id']
                )
                vir_genes.append(dc_asdict(gene))
            results[virulence_category] = vir_genes
    return results


@click.command()
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
def cli(sample_id, run_metadata, quast, process_metadata, kraken, mlst, cgmlst, virulence, resistance, output):
    """Combine pipeline results into a standardized json output file."""
    # base results
    results = {
        'sample_id': sample_id,
        'meta': {'run': json.load(run_metadata)},
        'qc': {},
    }
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