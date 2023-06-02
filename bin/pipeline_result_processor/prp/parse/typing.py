"""Parsers for various typing tools."""

import csv
import json
import logging

from ..models.sample import MethodIndex
from ..models.typing import TypingMethod, TypingResultCgMlst, TypingResultMlst
from ..models.typing import TypingSoftware as Software

LOG = logging.getLogger(__name__)


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
