"""Parse output of various phenotype predicting tools."""

import csv
import json
import logging
import pandas as pd
from typing import Any, Dict, Tuple

from ..models.metadata import SoupVersion, SoupVersions
from ..models.phenotype import (
    PhenotypeResult,
    ElementType,
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
        if (limit_to_phenotypes is not None and phenotype["key"] not in limit_to_phenotypes):
            continue

        if "resistant" in phenotype.keys():
            if phenotype["resistant"]:
                resistant.add(phenotype["resistance"])
            else:
                susceptible.add(phenotype["resistance"])
    return {"susceptible": list(susceptible), "resistant": list(resistant)}


def _parse_resfinder_amr_genes(resfinder_result, limit_to_phenotypes=None) -> Tuple[ResistanceGene, ...]:
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


def _parse_resfinder_amr_variants(resfinder_result, limit_to_phenotypes=None) -> Tuple[ResistanceVariant, ...]:
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

def _parse_amrfinder_amr_results(pred, element_type: str):
    """Parse amrfinder prediction results from amrfinderplus."""
    with open(pred, "rb") as tsvfile:
        hits = pd.read_csv(tsvfile, delimiter="\t")
        hits = hits.rename(columns={"Contig id": "contig_id", "Gene symbol": "gene_symbol", "Sequence name": "sequence_name", "Element type": "element_type", 
                                    "Element subtype": "element_subtype", "Target length": "target_length", "Reference sequence length": "ref_seq_len", 
                                    "% Coverage of reference sequence": "ref_seq_cov", "% Identity to reference sequence": "ref_seq_identity", 
                                    "Alignment length": "align_len", "Accession of closest sequence": "close_seq_accn", "Name of closest sequence": "close_seq_name"})
        hits = hits.drop(columns=["Protein identifier", "HMM id", "HMM description"])
        if element_type == "AMR":
            res = hits[hits["element_type"] == "AMR"].to_json(orient="records")
        elif element_type == "ENV":
            res = hits[(hits["element_subtype"] == "HEAT")].to_json(orient="records")
        elif element_type == "CHEM":
            res = hits[(hits["element_subtype"] == "ACID") & (hits["element_subtype"] == "BIOCIDE") & (hits["element_subtype"] == "METAL")].to_json(orient="records")
    return res


def parse_resfinder_amr_pred(prediction: Dict[str, Any], resistance_category) -> Tuple[SoupVersions, PhenotypeResult]:
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
        ElementType.CHEM: [
            "formaldehyde",
            "benzylkonium chloride",
            "ethidium bromide",
            "chlorhexidine",
            "cetylpyridinium chloride",
            "hydrogen peroxide",
        ],
        ElementType.ENV: ["temperature"],
    }
    categories[ElementType.AMR] = list(
        {k for k in prediction["phenotypes"].keys()}
        - set(categories[ElementType.CHEM] + categories[ElementType.ENV])
    )

    # parse resistance
    resistance = PhenotypeResult(
        phenotypes = _get_resfinder_amr_sr_profie(prediction, categories[resistance_category]),
        genes = _parse_resfinder_amr_genes(prediction, categories[resistance_category]),
        mutations = _parse_resfinder_amr_variants(prediction, categories[resistance_category]),
    )
    return MethodIndex(type=resistance_category, result=resistance)


def parse_amrfinder_amr_pred(file, element_type: str):
    """Parse amrfinder resistance prediction results."""
    LOG.info("Parsing amrfinder amr prediction")
    if file:
        results = _parse_amrfinder_amr_results(file, element_type)
    else:
        results = _default_resistance()
    return results


def _parse_virulencefinder_vir_results(pred: str) -> PhenotypeResult:
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

def _parse_amrfinder_vir_results(pred: str):
    """Parse amrfinder prediction results from amrfinderplus."""
    hits = pd.read_csv(pred, delimiter="\t")
    hits = hits.rename(columns={"Contig id": "contig_id", "Gene symbol": "gene_symbol", "Sequence name": "sequence_name", "Element type": "element_type", 
                                "Element subtype": "element_subtype", "Target length": "target_length", "Reference sequence length": "ref_seq_len", 
                                "% Coverage of reference sequence": "ref_seq_cov", "% Identity to reference sequence": "ref_seq_identity", 
                                "Alignment length": "align_len", "Accession of closest sequence": "close_seq_accn", "Name of closest sequence": "close_seq_name"})
    hits = hits.drop(columns=["Protein identifier", "HMM id", "HMM description"])
    virulence = hits[hits["element_type"] == "VIRULENCE"].to_json(orient="records")
    return virulence

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


def parse_virulencefinder_vir_pred(file: str) -> PhenotypeResult:
    """Parse virulencefinder virulence prediction results."""
    LOG.info("Parsing virulencefinder virulence prediction")
    pred = json.load(file)
    if "not virulencefinder" in pred:
        results: PhenotypeResult = _parse_virulencefinder_vir_results(pred)
    else:
        results: PhenotypeResult = _default_virulence()
    return MethodIndex(type=ElementType.VIR, result=results)

def parse_amrfinder_vir_pred(file: str):
    """Parse amrfinder virulence prediction results."""
    LOG.info("Parsing amrfinder virulence prediction")
    results = _parse_amrfinder_vir_results(file)
    return results
