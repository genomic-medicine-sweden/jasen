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
            return int(allele.split("-")[1])
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


def _index_to_record(idx_dict):
    rec_dict = []
    for lineage in idx_dict:
        variant = list(idx_dict[lineage].keys())[0]
        lin_array = {
            "lin": lineage,
            "family": None,
            "spoligotype": None,
            "rd": None,
            "frac": None,
            "variant": variant,
            "coverage": idx_dict[lineage][variant]["info"]["coverage"]["alternate"],
        }
        rec_dict.append(lin_array)
    return rec_dict


def parse_tbprofiler_lineage_results(pred_res: dict, method) -> TypingResultLineage:
    """Parse tbprofiler results for lineage object."""
    LOG.info("Parsing lineage results")
    #lineages = pd.read_json(json.dumps(pred_res["lineage"]), orient="records").to_json(orient="index", index="lin")
    #print(lineages)
    result_obj = TypingResultLineage(
        main_lin=pred_res["main_lin"],
        sublin=pred_res["sublin"],
        lineages=pred_res["lineage"],
    )
    return MethodIndex(type=method, software=Software.TBPROFILER, result=result_obj)

def parse_mykrobe_lineage_results(pred_res: dict, method) -> TypingResultLineage:
    """Parse mykrobe results for lineage object."""
    LOG.info("Parsing lineage results")
    try:
        genotypes = list(list(pred_res["phylogenetics"]["lineage"]["calls_summary"].values())[0]["genotypes"].keys())
        lineage_calls = pred_res["phylogenetics"]["lineage"]["calls"]
        sublin = list(lineage_calls.keys())[0]
        result_obj = TypingResultLineage(
            main_lin=genotypes[0],
            sublin=sublin,
            lineages=_index_to_record(lineage_calls[sublin]),
        )
    except KeyError:
        genotypes = list(pred_res["phylogenetics"]["lineage"].keys())
        result_obj = TypingResultLineage(
            main_lin=genotypes[0],
            sublin=genotypes[0],
            lineages=[{
                "lin": genotypes[0],
                "family": None,
                "spoligotype": None,
                "rd": None,
                "frac": None,
                "variant": None,
                "coverage": pred_res["phylogenetics"]["lineage"][genotypes[0]],
            }],
        )
    return MethodIndex(type=method, software=Software.MYKROBE, result=result_obj)
