#! /usr/bin/env python
"""Combine output files from pipeline into a standardized json output."""

import csv
import json
import logging
from collections import defaultdict
from typing import Any, Dict, Tuple
from unittest import result

import click
import pandas as pd

from .models.metadata import SoupVersion, SoupVersions
from .models.phenotype import (
    PhenotypeResult,
    ResistanceGene,
    ResistanceVariant,
    VirulenceGene,
    PhenotypeType,
)
from .models.sample import AssemblyQc, MethodIndex
from .models.typing import TypingResultCgMlst, TypingResultMlst, TypingMethod

LOG = logging.getLogger(__name__)


SPP_MIN_READ_FRAC = 0.001


def parse_qust_results(file: str) -> AssemblyQc:
    """Parse quast file and extract relevant metrics.

    Args:
        sep (str): seperator

    Returns:
        AssemblyQc: list of key-value pairs
    """
    LOG.info(f"Parsing tsv file: {file}")
    creader = csv.reader(file, delimiter="\t")
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
        duplication_ratio=raw[0]["Duplication ratio"],
    )


def parse_mlst_results(path: str) -> TypingResultMlst:
    """Parse mlst results from mlst to json object."""
    LOG.info("Parsing mlst results")
    result = json.load(path)[0]
    result_obj = TypingResultMlst(
        scheme=result["scheme"],
        sequence_type=None
        if result["sequence_type"] == "-"
        else result["sequence_type"],
        alleles={
            gene: None if allele == "-" else allele
            for gene, allele in result["alleles"].items()
        },
    )
    return MethodIndex(type=TypingMethod.mlst, result=result_obj)


def parse_cgmlst_results(
    file: str, include_novel_alleles: bool = True
) -> TypingResultCgMlst:
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
    ERRORS = ["LNF", "PLOT5", "NIPH", "NIPHEM", "ALM", "ASM"]

    def replace_errors(allele):
        """Replace errors and novel alleles with nulls if they are not to be inlcuded."""
        if any(
            [allele in ERRORS, allele.startswith("INF") and not include_novel_alleles]
        ):
            return None
        elif allele.startswith("INF") and include_novel_alleles:
            return int(allele.split("-")[1])
        return int(allele)

    msg = "Parsing cgmslt results, "
    LOG.info(
        msg + "not" if not include_novel_alleles else "" + "including novel alleles"
    )
    creader = csv.reader(file, delimiter="\t")
    _, *allele_names = (colname.rstrip(".fasta") for colname in next(creader))
    # parse alleles
    _, *alleles = next(creader)
    corrected_alleles = (replace_errors(a) for a in alleles)
    results = TypingResultCgMlst(
        n_novel=sum(1 for a in alleles if a.startswith("INF")),
        n_missing=sum(1 for a in alleles if a in ERRORS),
        alleles=dict(zip(allele_names, corrected_alleles)),
    )
    return MethodIndex(type=TypingMethod.cgmlst, result=results)


def _get_resfinder_amr_sr_profie(resfinder_result, limit_to_phenotypes=None):
    """Get resfinder susceptibility/resistance profile."""
    susceptible = set()
    resistant = set()
    for phenotype in resfinder_result["phenotypes"].values():
        # skip phenotype if its not part of the desired category
        if (
            limit_to_phenotypes is not None
            and phenotype["key"] not in limit_to_phenotypes
        ):
            continue

        if phenotype["resistant"]:
            resistant.add(phenotype["resistance"])
        else:
            susceptible.add(phenotype["resistance"])
    return {"susceptible": list(susceptible), "resistant": list(resistant)}


def _parse_resfinder_amr_genes(
    resfinder_result, limit_to_phenotypes=None
) -> Tuple[ResistanceGene, ...]:
    """Get resistance genes from resfinder result."""
    results = []
    for info in resfinder_result["genes"].values():
        # Get only acquired resistance genes
        if not info["ref_database"].startswith("Res"):
            continue

        # Get only genes of desired phenotype
        if limit_to_phenotypes is not None:
            intersect = set(info["phenotypes"]) & set(limit_to_phenotypes)
            if len(intersect) == 0:
                continue

        # store results
        gene = ResistanceGene(
            name=info["name"],
            accession=info["ref_acc"],
            depth=info["depth"],
            identity=info["identity"],
            coverage=info["coverage"],
            ref_start_pos=info["ref_start_pos"],
            ref_end_pos=info["ref_end_pos"],
            ref_gene_length=info["ref_gene_lenght"],
            alignment_length=info["alignment_length"],
            phenotypes=info["phenotypes"],
            ref_database=info["ref_database"],
            ref_id=info["ref_id"],
        )
        results.append(gene)
    return results


def _parse_resfinder_amr_variants(
    resfinder_result, limit_to_phenotypes=None
) -> Tuple[ResistanceVariant, ...]:
    """Get resistance genes from resfinder result."""
    results = []
    for info in resfinder_result["seq_variations"].values():
        # Get only variants from desired phenotypes
        if limit_to_phenotypes is not None:
            intersect = set(info["phenotypes"]) & set(limit_to_phenotypes)
            if len(intersect) == 0:
                continue
        # get gene depth
        info["depth"] = resfinder_result["genes"][info["genes"][0]]["depth"]
        # translate variation type bools into classifier
        if info["substitution"]:
            var_type = "substitution"
        elif info["insertion"]:
            var_type = "insertion"
        elif info["deletion"]:
            var_type = "deletion"
        else:
            raise ValueError("Output has no known mutation type")
        variant = ResistanceVariant(
            variant_type=var_type,
            genes=info["genes"],
            phenotypes=info["phenotypes"],
            position=info["ref_start_pos"],
            ref_codon=info["ref_codon"],
            alt_codon=info["var_codon"],
            depth=info["depth"],
            ref_database=info["ref_database"],
            ref_id=info["ref_id"],
        )
        results.append(variant)
    return results


def parse_resistance_pred(
    prediction: Dict[str, Any], resistance_category
) -> Tuple[SoupVersions, PhenotypeResult]:
    """Parse resfinder resistance prediction results."""
    # resfinder missclassifies resistance the param amr_category by setting all to amr
    LOG.info("Parsing resistance prediction")
    meta = [
        SoupVersion(
            **{
                "name": metrics["database_name"],
                "version": metrics["database_version"],
                "type": metrics["type"],
            }
        )
        for metrics in prediction["databases"].values()
    ]
    # parse resistance based on the category
    categories = {
        PhenotypeType.chem: [
            "formaldehyde",
            "benzylkonium chloride",
            "ethidium bromide",
            "chlorhexidine",
            "cetylpyridinium chloride",
            "hydrogen peroxide",
        ],
        PhenotypeType.env: ["temperature"],
    }
    categories[PhenotypeType.amr] = list(
        {k for k in prediction["phenotypes"].keys()}
        - set(categories[PhenotypeType.chem] + categories[PhenotypeType.env])
    )

    # parse resistance
    resistance = PhenotypeResult(
        phenotypes=_get_resfinder_amr_sr_profie(
            prediction, categories[resistance_category]
        ),
        genes=_parse_resfinder_amr_genes(prediction, categories[resistance_category]),
        mutations=_parse_resfinder_amr_variants(
            prediction, categories[resistance_category]
        ),
    )
    return MethodIndex(type=resistance_category, result=resistance)


def _parse_virulence_finder_results(pred: str) -> PhenotypeResult:
    """Parse virulence prediction results from ARIBA."""
    results = {}
    # parse virulence finder results
    species = [k for k in pred["virulencefinder"]["results"].keys()]
    for key, genes in pred["virulencefinder"]["results"][species[0]].items():
        virulence_category = key.split("_")[1]
        vir_genes = []
        for gn in genes.values():
            start_pos, end_pos = map(int, gn["position_in_ref"].split(".."))
            gene = VirulenceGene(
                name=gn["virulence_gene"],
                virulence_category=virulence_category,
                accession=gn["accession"],
                depth=None,
                identity=gn["identity"],
                coverage=gn["coverage"],
                ref_start_pos=start_pos,
                ref_end_pos=end_pos,
                ref_gene_length=gn["template_length"],
                alignment_length=gn["HSP_length"],
                ref_database="virulenceFinder",
                ref_id=gn["hit_id"],
            )
            vir_genes.append(gene)
        results[virulence_category] = vir_genes
    return PhenotypeResult(results)


def _parse_ariba_results(pred: str) -> PhenotypeResult:
    """Parse virulence prediction results from ARIBA."""
    absent_genes = []
    present_genes = []
    for gene, metrics in pred.items():
        if metrics["present"] == 0:
            absent_genes.append(gene)
            continue
        best_hit = max(
            [k for k in metrics.keys() if not k == "present"],
            key=lambda x: int(x.split("_")[-1]),
        )
        hit = metrics[best_hit]
        vir_gene = VirulenceGene(
            name=gene,
            virulence_category="",
            accession="",
            depth=None,
            identity=hit["id"],
            coverage=int(hit["match_len"]) / int(hit["ref_len"]),
            ref_start_pos=0,
            ref_end_pos=0,
            ref_gene_length=int(hit["ref_len"]),
            alignment_length=int(hit["match_len"]),
            ref_database="ariba",
            ref_id=best_hit,
        )
        present_genes.append(vir_gene)
    return PhenotypeResult(phenotypes=[], genes=present_genes, mutations=[])


def parse_virulence_pred(file: str) -> PhenotypeResult:
    """Parse virulence prediction results."""
    LOG.info("Parsing virulence prediction")
    pred = json.load(file)
    if "virulencefinder" in pred:
        results: PhenotypeResult = _parse_virulence_finder_results(pred)
    elif "ariba" in pred:
        results: PhenotypeResult = _parse_ariba_results(pred)
    else:
        results: PhenotypeResult = _parse_ariba_results(pred)

    return MethodIndex(type=PhenotypeType.vir, result=results)


def parse_species_pred(file: str):
    """ "parse_species_pred""Parse species prediciton result"""
    specie_pred: pd.DataFrame = pd.read_csv(file, sep="\t").sort_values(
        "fraction_total_reads", ascending=False
    )
    # limit the number of predicted species
    specie_pred = specie_pred[specie_pred["fraction_total_reads"] > SPP_MIN_READ_FRAC]
    import pdb

    pdb.set_trace()
    return specie_pred.to_dict(orient="records")
