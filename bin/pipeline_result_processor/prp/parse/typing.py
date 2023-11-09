"""Parsers for various typing tools."""

import csv
import json
import logging
import pandas as pd

from ..models.sample import MethodIndex
from ..models.typing import TypingMethod, TypingResultCgMlst, TypingResultMlst, TypingResultMlst, TypingResultLineage
from ..models.typing import TypingSoftware as Software

LOG = logging.getLogger(__name__)

def _process_allele_call(allele):
    if allele.isdigit():
        return int(allele)
    elif ',' in allele:
        return allele.split(',')
    elif '?' in allele:
        return "partial"
    elif '~' in allele:
        return "novel"
    elif allele == '-':
        return
    raise ValueError(f"MLST allele {allele} not expected format")


def parse_mlst_results(path: str) -> TypingResultMlst:
    """Parse mlst results from mlst to json object."""
    LOG.info("Parsing mlst results")
    result = json.load(path)[0]
    result_obj = TypingResultMlst(
        scheme=result["scheme"],
        sequence_type=None
        if result["sequence_type"] == "-"
        else result["sequence_type"],
        alleles={gene: _process_allele_call(allele) for gene, allele in result["alleles"].items()},
    )
    return MethodIndex(
        type=TypingMethod.MLST, software=Software.MLST, result=result_obj
    )


def parse_cgmlst_results(
    file: str, include_novel_alleles: bool = True, correct_alleles: bool = False
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
    ERRORS = ("LNF", "PLOT3", "PLOT5", "NIPH", "NIPHEM", "ALM", "ASM")

    def replace_errors(allele):
        """Replace errors and novel alleles with nulls if they are not to be inlcuded."""
        if any(
            [
                correct_alleles and allele in ERRORS,
                correct_alleles
                and allele.startswith("INF")
                and not include_novel_alleles,
            ]
        ):
            return None
        elif allele.startswith("INF") and include_novel_alleles:
            try:
                allele = int(allele.split("-")[1])
            except ValueError:
                allele = str(allele.split("-")[1])
        try:
            allele = int(allele)
        except ValueError:
            allele = str(allele)
        return allele

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
    return MethodIndex(
        type=TypingMethod.CGMLST, software=Software.CHEWBBACA, result=results
    )


def _create_lineage_array(lineage=None, family=None, spoligotype=None, rd=None, frac=None, variant=None, coverage=None):
    """Create lineage array for mykrobe info"""
    return {
        "lin": lineage,
        "family": family,
        "spoligotype": spoligotype,
        "rd": rd,
        "frac": frac,
        "variant": variant,
        "coverage": coverage,
    }


def _get_lineage_info(lineage_dict):
    """Create a list of arrays by parsing mykrobe lineage output"""
    lineages = []
    if "calls_summary" in lineage_dict:
        lineage_calls = lineage_dict["calls"]
        sublin = list(lineage_calls.keys())[0]
        lineage_info = lineage_calls[sublin]
        for lineage in lineage_info:
            genotypes = list(list(lineage_dict["calls_summary"].values())[0]["genotypes"].keys())
            main_lin = genotypes[0]
            try:
                variant = list(lineage_info[lineage].keys())[0]
            except AttributeError:
                variant = None
            try:
                coverage = lineage_info[lineage][variant]["info"]["coverage"]["alternate"]
            except (KeyError, TypeError):
                coverage = None
            lin_array = _create_lineage_array(lineage=lineage, variant=variant, coverage=coverage)
            lineages.append(lin_array)
    else:
        genotypes = list(lineage_dict.keys())
        main_lin, sublin = genotypes[0], genotypes[0]
        lin_array = _create_lineage_array(lineage=genotypes[0], coverage=lineage_dict[genotypes[0]])
        lineages.append(lin_array)
    return main_lin, sublin, lineages


def parse_tbprofiler_lineage_results(pred_res: dict, method) -> TypingResultLineage:
    """Parse tbprofiler results for lineage object."""
    LOG.info("Parsing lineage results")
    result_obj = TypingResultLineage(
        main_lin=pred_res["main_lin"],
        sublin=pred_res["sublin"],
        lineages=pred_res["lineage"],
    )
    return MethodIndex(type=method, software=Software.TBPROFILER, result=result_obj)


def parse_mykrobe_lineage_results(pred_res: dict, method) -> TypingResultLineage:
    """Parse mykrobe results for lineage object."""
    LOG.info("Parsing lineage results")
    main_lin, sublin, lineages = _get_lineage_info(pred_res["phylogenetics"]["lineage"])
    result_obj = TypingResultLineage(
        main_lin=main_lin,
        sublin=sublin,
        lineages=lineages,
    )
    return MethodIndex(type=method, software=Software.MYKROBE, result=result_obj)
