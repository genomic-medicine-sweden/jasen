"""Parse output of QC tools."""

import csv
import logging

from ..models.qc import QcMethodIndex, QcSoftware, QuastQcResult
from click.types import File

LOG = logging.getLogger(__name__)


def parse_quast_results(file: File) -> QcMethodIndex:
    """Parse quast file and extract relevant metrics.

    Args:
        sep (str): seperator

    Returns:
        AssemblyQc: list of key-value pairs
    """
    LOG.info(f"Parsing tsv file: {file.name}")
    creader = csv.reader(file, delimiter="\t")
    header = next(creader)
    raw = [dict(zip(header, row)) for row in creader]
    qc_res = QuastQcResult(
        total_length=int(raw[0]["Total length"]),
        reference_length=raw[0]["Reference length"],
        largest_contig=raw[0]["Largest contig"],
        n_contigs=raw[0]["# contigs"],
        n50=raw[0]["N50"],
        assembly_gc=raw[0]["GC (%)"],
        reference_gc=raw[0]["Reference GC (%)"],
        duplication_ratio=raw[0]["Duplication ratio"],
    )
    return QcMethodIndex(software=QcSoftware.QUAST, result=qc_res)
