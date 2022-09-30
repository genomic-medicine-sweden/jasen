"""Parse output of various phenotype predicting tools."""

import csv
import json
import logging
from typing import Any, Dict, Tuple

from ..models.metadata import SoupVersion, SoupVersions
from ..models.phenotype import (
    PhenotypeResult,
    PhenotypeType,
    ResistanceGene,
    ResistanceVariant,
    VirulenceGene,
)
from ..models.sample import MethodIndex

LOG = logging.getLogger(__name__)


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

        if "resistant" in phenotype.keys():
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

    if not "genes" in resfinder_result.keys():
        results  = _default_resistance().genes
        return results

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
    igenes = []
    for info in resfinder_result["seq_variations"].values():
        # Get only variants from desired phenotypes
        if limit_to_phenotypes is not None:
            intersect = set(info["phenotypes"]) & set(limit_to_phenotypes)
            if len(intersect) == 0:
                continue
        # get gene depth
        if "genes" in resfinder_result.keys():
            info["depth"] = resfinder_result["genes"][info["genes"][0]]["depth"]
        else:
            info["depth"] = 0
        # translate variation type bools into classifier
        if info["substitution"]:
            var_type = "substitution"
        elif info["insertion"]:
            var_type = "insertion"
        elif info["deletion"]:
            var_type = "deletion"
        else:
            raise ValueError("Output has no known mutation type")
        if not "genes" in info.keys():
            #igenes = _default_resistance().genes
            igenes = [""]
        variant = ResistanceVariant(
            variant_type=var_type,
            genes=igenes,
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
        PhenotypeType.CHEM: [
            "formaldehyde",
            "benzylkonium chloride",
            "ethidium bromide",
            "chlorhexidine",
            "cetylpyridinium chloride",
            "hydrogen peroxide",
        ],
        PhenotypeType.ENV: ["temperature"],
    }
    categories[PhenotypeType.AMR] = list(
        {k for k in prediction["phenotypes"].keys()}
        - set(categories[PhenotypeType.CHEM] + categories[PhenotypeType.ENV])
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
        if not genes == "No hit found":
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


def _default_virulence() -> PhenotypeResult:
    gene = VirulenceGene(
        name="none",
        virulence_category="",
        accession="",
        depth=None,
        identity=0,
        coverage=0,
        ref_start_pos=0,
        ref_end_pos=0,
        ref_gene_length=0,
        alignment_length=0,
        ref_database="",
        ref_id=0,
    )
    genes = list()
    genes.append(gene)
    return PhenotypeResult(phenotypes=[], genes=genes, mutations=[])


def _default_resistance() -> PhenotypeResult:
    gene = ResistanceGene(
        name="none",
        virulence_category="",
        accession="",
        depth=None,
        identity=0,
        coverage=0,
        ref_start_pos=0,
        ref_end_pos=0,
        ref_gene_length=0,
        alignment_length=0,
        ref_database="",
        phenotypes=[],
        ref_id=0,
    )
    genes = list()
    genes.append(gene)
    return PhenotypeResult(phenotypes=[], genes=genes, mutations=[])


def parse_virulence_pred(file: str) -> PhenotypeResult:
    """Parse virulence prediction results."""
    LOG.info("Parsing virulence prediction")
    pred = json.load(file)
    if "not virulencefinder" in pred:
        results: PhenotypeResult = _parse_virulence_finder_results(pred)
    elif "ariba" in pred:
        results: PhenotypeResult = _parse_ariba_results(pred)
    else:
        results: PhenotypeResult = _default_virulence()
    return MethodIndex(type=PhenotypeType.VIR, result=results)
